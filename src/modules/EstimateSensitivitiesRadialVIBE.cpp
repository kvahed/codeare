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

	return codeare::OK;
}

codeare::error_code EstimateSensitivitiesRadialVIBE::Prepare () {
	Matrix<cxfl> sensitivities;
	Matrix<float> kspace, weights;
	Add<cxfl>("sensitivities", sensitivities);
	Add<float>("kspace", kspace);
	Add<float>("weights", weights);

	return codeare::OK;
}

codeare::error_code EstimateSensitivitiesRadialVIBE::Process () {
	// Matrices
	Matrix<cxfl>& meas = Get<cxfl>("meas");
	meas = squeeze(1.e8*meas);
	Matrix<cxfl>& sensitivities = Get<cxfl>("sensitivities");
	Matrix<float>& kspace = Get<float>("kspace");
	Matrix<float>& weights = Get<float>("weights");
	Matrix<cxfl> orig = meas;
	size_t nk = size(meas,0), nc = size(meas,1), nv = size(meas,2), nz = size(meas,3);

	// slices, samples, coils, projections [3 0 1 2]
	Vector<size_t> perm(4);
	perm[0]=3; perm[1]=0; perm[2]=1; perm[3]=2;
	meas = permute(meas,perm);

	// FFT slice dimension
	// permute dims for NuFFT
	// samples, projections, slices, coils [1 3 0 2]
	perm[0]=1; perm[1]=3; perm[2]=0; perm[3]=2;
	meas = permute(ifft(meas), perm);

	// Remove top and bottom slices
	meas = meas(CR(),CR(),CR(2,nz-3),CR()); nz -= 4;
	meas = resize(meas,nv*nk,nz,nc);

	// Sort for channel NuFFT: samples, slices,
	kspace = zeros<float>(2,nk,nv);
	kspace(R(0),R(),R(0)) = linspace<float>(-.5, .5, nk);
	cxfl rot = std::polar<float>(1.0f,NYUGA_RAD);
	for (size_t j = 1; j < nv; ++j)
		for (size_t i = 0; i < nk; ++i) {
			cxfl tmp = rot * cxfl(kspace(0,i,j-1),kspace(1,i,j-1));
			kspace(0,i,j) = std::real(tmp);
			kspace(1,i,j) = std::imag(tmp);
		}
	kspace = resize(kspace,2,nv*nk);
	weights = ones<float>(nv*nk,1);

	Vector<size_t> sens_dims = _image_space_dims;
	sens_dims.PushBack(nc);
	sensitivities = Matrix<cxfl>(sens_dims);
	Matrix<float> density_comp = zeros<float>(nk,1);

	// FT operators
	Vector<NFFT<cxfl> > FTOP;
	Params p;
	p["nk"]           = 4096;
	p["imsz"]         = _image_space_dims;
	p["3rd_dim_cart"] = _cart_3rd_dim;
	FTOP.PushBack(NFFT<cxfl>(p));
	FTOP[0].KSpace(kspace(CR(),CR(0,4095)));
	FTOP[0].Weights(weights(CR(0,4095)));
	std::cout << FTOP[0] << std::endl;
	sensitivities (R(),R(),R(),R(0)) = FTOP[0] ->* meas(CR(0,4095),CR(),CR(0));
	Add<cxfl>("orig", orig);
	Add<float>("weights", weights);

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
