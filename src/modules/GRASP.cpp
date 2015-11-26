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

#include "GRASP.hpp"
#include "Toolbox.hpp"
#ifdef HAVE_NFFT3
	#include "NCSENSE.hpp"
#endif
#include "CS_XSENSE.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

using namespace RRStrategy;

codeare::error_code GRASP::Init () {

	printf ("Intialising GRASP ...\n");

    std::cout << *this << std::endl;

	for (size_t i = 0; i < 3; i++)
		m_N[i] = 1;



	m_initialised = true;

	printf ("... done.\n\n");

	return codeare::OK;

}

codeare::error_code GRASP::Prepare () {
	codeare::error_code error = codeare::OK;
	return error;
}

codeare::error_code GRASP::Process () {

    Matrix<cxfl>& data = Get<cxfl>("meas");
    Matrix<float>& kspace = Get<float>("kspace");
    Matrix<float>& weights = Get<float>("weights");

    size_t nline, nv, nt, nx, nz, nc;
    float ta = wspace.PGet<float>("TA");
    Vector<size_t> n = size(data);
    nx = n[0]; nv = n[1]; nz = n[2]; nc = n[3];

    // Timing
	nt = ceil(ta/_tf);
	nv = floor (nv/nt);

	std::cout << "  Data acquired within " << ta << std::endl;
	std::cout << "  Sorting data in " << nt << " contrast gates with "	<< nv << " views each." <<std::endl;

	// Clip over the time data, kspace and signal
	data = data(CR(),CR(0,nt*nv-1),CR(),CR());
	kspace = kspace(CR(),CR(),CR(0,nt*nv-1));
	nv = size(data,1);

	std::cout << "  Incoming    " << std::endl;
	std::cout << "    data:     " << size(data) << std::endl;
	std::cout << "    kspace:   " << size(kspace) << std::endl;

	// Contrasts sorting
	data = resize(data,nx,nv/nt,nt,nz,nc);
	kspace = resize(kspace,size(kspace,0),nx,nv/nt,nt);

	std::cout << "  Reshaped:   " << std::endl;
	std::cout << "    data:     " << size(data) << std::endl;
	std::cout << "    kspace:   " << size(kspace) << std::endl;

	Vector<size_t> order(4); order[0]=0; order[1]=2; order[2]=3; order[3]=1;
	data = permute (data,order);
	kspace = resize(kspace,size(kspace,0),nx*nv,nt);
	weights = zeros<float>(nx*nv,1);
    weights (R( 0,nx/2-1),0) = linspace<float>(1.,1./nx,nx/2);
    weights (R(nx/2,nx-1),0) = linspace<float>(1./nx,1.,nx/2);
    for (auto i = weights.Begin()+nx; i < weights.End(); i += nx)
        std::copy (weights.Begin(), weights.Begin()+nx, i);

	std::cout << "  Reshaped and permuted:    " << std::endl;
	std::cout << "    data:     " << size(data) << std::endl;
	std::cout << "    kspace:   " << size(kspace) << std::endl;
	std::cout << "    weights:  " << size(weights) << std::endl;

	// FT operator
	Matrix<cxfl>& sensitivities = Get<cxfl>("sensitivities");
	ft_params["sensitivities"] = sensitivities;
	Vector<size_t> image_size = size(sensitivities); image_size.pop_back();
	ft_params["imsz"] = image_size;
	ft_params["dim4"] = (int) nt;
	ft_params["nk"]   = size(data,0);
    CS_XSENSE<cxfl> ft(ft_params);
    ft.KSpace(kspace);
	ft.Weights (weights);
	std::cout << ft << std::endl;

	// 1-3 NuFFT 4 TV CS
	Matrix<cxfl> im_xd = ft ->* data;

    Add ("im_xd", im_xd);

    return codeare::OK;

}


GRASP::GRASP() : m_wm(0), m_csiter(0), m_wf(0), m_dim(0), m_verbose(0), m_noise(0.),
		m_test_case(0), _tf(15.0) {}


GRASP::~GRASP() {}


codeare::error_code GRASP::Finalise() {return codeare::OK;}


// the class facxflories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new GRASP;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}


