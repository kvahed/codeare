/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but 
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301  USA
 */

#include "Algos.hpp"
#include "EstimateSensitivitiesRadialVIBE.hpp"
#include "DFT.hpp"
#include "NFFT.hpp"
#include "Creators.hpp"


static const float NYUGA_RAD = PI*111.246117975f/180.0f;

using namespace RRStrategy;

codeare::error_code EstimateSensitivitiesRadialVIBE::Init () {
	_image_space_dims = GetList<size_t>("image_space_dims");
	_cart_3rd_dim     = GetAttr<bool>("cart_3rd_dim");
	Matrix<cxfl> sensitivities;
	Add ("sensitivities", sensitivities);
	return codeare::OK;
}

codeare::error_code EstimateSensitivitiesRadialVIBE::Prepare () {
	Matrix<float> kspace, weights, sos;
	Add<float>("kspace", kspace);
	Add<float>("weights", weights);
	Add<float>("sos", sos);
	return codeare::OK;
}

inline void EstimateSensitivitiesRadialVIBE::FormGARadialKSpace (
		const size_t& nk, const size_t& nv) const {
	Matrix<float>& kspace = Get<float>("kspace");
	Matrix<float>& weights = Get<float>("weights");
	kspace = zeros<float>(2,nk,nv);
	kspace(R(0),R(),R(0)) = linspace<float>(-.45, .45, nk);
	cxfl rot = std::polar<float>(1.0f,NYUGA_RAD);
	for (size_t j = 1; j < nv; ++j)
		for (size_t i = 0; i < nk; ++i) {
			cxfl tmp = rot * cxfl(kspace(0,i,j-1),kspace(1,i,j-1));
			kspace(0,i,j) = std::real(tmp);
			kspace(1,i,j) = std::imag(tmp);
		}
	kspace = resize(kspace,2,nv*nk);
	weights = ones<float>(nv*nk,1);
}


codeare::error_code EstimateSensitivitiesRadialVIBE::Process () {
	// Matrices
	Matrix<cxfl>& meas = Get<cxfl>("meas");
	meas = squeeze(2.e6*meas);
	Matrix<cxfl>& sensitivities = Get<cxfl>("sensitivities");
	Matrix<float>& kspace = Get<float>("kspace");
	Matrix<float>& weights = Get<float>("weights");
	Matrix<float>& sos = Get<float>("sos");
	std::cout << "  Incoming: " << size(meas) << std::endl;
	size_t nk = size(meas,0), nv = size(meas,1), nz = size(meas,2), nc = size(meas,3);

	// Permute for slice FFT [2 0 1 3] & FFT
	Vector<size_t> perm(4);
	perm[0]=2; perm[1]=0; perm[2]=1; perm[3]=3;
	meas = permute(meas,perm);
	std::cout << "  Permuted for slice FFT: " << size(meas) << std::endl;

	// Slice FFT
	std::cout << "  Performing slice FFT ..." << std::endl;
	meas = ifft(meas,0,true);

	// Removing top and bottom slices
    size_t pad = (nz - _image_space_dims[2])/2;
    meas = meas(CR(7,nz-8),CR(),CR(),CR()); nz = size(meas,0);
	std::cout << "  Slice direction reduced: " << size(meas) << std::endl;

	// Permute for global coil nufft [1 2 0 3]
	perm[0]=1; perm[1]=2; perm[2]=0; perm[3]=3;
	meas = permute(meas, perm);
	std::cout << "  Permuted for channel NuFFTs: " << size(meas) << std::endl;
	meas = resize(meas,nv*nk,nz,nc);
	std::cout << "  Collapsed samples and views dimensions: " << size(meas) << std::endl;

	// Sort for channel NuFFT: samples, slices,
	std::cout << "  Building GA stack of star k-space trajectory ..." << std::endl;
	FormGARadialKSpace(nk, nv);

	// FT operators
	std::cout << "  Building NuFFT operator(s) ..." << std::endl;
	Vector<size_t> sens_dims = _image_space_dims;
	sens_dims.push_back(nc);
	sensitivities = Matrix<cxfl>(sens_dims);
	Matrix<float> density_comp = zeros<float>(nk,1);
	Vector<NFFT<cxfl> > FTOP;
	Params p;
	p["nk"] = nv*nk; p["imsz"] = _image_space_dims; p["3rd_dim_cart"] = _cart_3rd_dim;
	p["m"] = (size_t)1; p["alpha"] = 1.0f; p["epsilon"]=7.e-4f; p["maxit"]=(size_t)2;
    size_t threads = 1;
    for (size_t i = 0; i < threads; ++i) {
        FTOP.push_back(NFFT<cxfl>(p));
        FTOP[i].KSpace(kspace);
        FTOP[i].Weights(weights);
    }
	std::cout << FTOP[0] << std::endl;
    
	// Channel NuFFTs
	std::cout << "  NuFFTing ..." << std::endl;
#pragma omp parallel for default (shared) num_threads(threads)
    for (size_t i = 0; i < nc; ++i) 
        sensitivities (R(),R(),R(),R(i)) = FTOP[omp_get_thread_num()] ->* meas(CR(),CR(),CR(i));

    sos = sum(abs(sensitivities),3)+1.e-9;
    
    for (size_t i = 0; i < nc; ++i)
        sensitivities (R(),R(),R(),R(i)) /= sos;

    meas = resize(meas,nk,nv,nz,nc);
    kspace = resize(kspace,size(kspace,0),nk,nv);

	return codeare::OK;
}

codeare::error_code EstimateSensitivitiesRadialVIBE::Finalise () {
	return codeare::OK;
}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                  {
    return new EstimateSensitivitiesRadialVIBE;
}
extern "C" DLLEXPORT void           destroy (ReconStrategy* p)  {
	delete p;
}
