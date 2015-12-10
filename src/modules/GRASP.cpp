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

enum GRASP_EXCEPTION {
		OK
};


GRASP::GRASP() {
	ft_params["cgeps"]   = 1.0e-6f;
	ft_params["csiter"]  = 3;
	ft_params["cgconv"]  = 0.0e-3f;
	ft_params["lsiter"]  = 6;
	ft_params["verbose"] = false;
    ft_params["xfmw"]    = 0.0f;
	ft_params["m"]       = (size_t)1;
	ft_params["cgiter"]  = (size_t)0;
	ft_params["lambda"]  = 1.0e-9f;
	ft_params["tvw1"]    = 0.02f;
	ft_params["l1"]      = 1.0e-15f;
	ft_params["pnorm"]   = 1.0f;
	ft_params["verbose"] = 0;
	ft_params["lsa"]     = 0.01f; // Linesearch lower limit
	ft_params["nlopt"]   = 0;
	ft_params["lsb"]     = 0.6f; // Linesearch upper limit
	ft_params["alpha"]   = 1.0f; // NuFFT oversampling
	ft_params["fteps"]   = 7.e-4f; // NuFFT gridding convergence
	ft_params["parallel_linesearch"]   = true; // Line search parallelism
	ft_params["ft"]      = 3; // NC SENSE
	ft_params["ftiter"]  = 1; // NuFFT gridding iterations
	_tf                  = 15;
	ft_params["nliter"]  = 6;
}

codeare::error_code GRASP::Init () {

	printf ("Intialising GRASP ...\n");
	try {
		ft_params["csiter"] = GetAttr<int>("csiter");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["nliter"] = GetAttr<int>("nliter");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["lsiter"] = GetAttr<int>("lsiter");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["cgiter"] = GetAttr<size_t>("cgmaxit");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["tvw1"] = GetAttr<float>("tvw");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["l1"] = GetAttr<float>("l1");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["pnorm"] = GetAttr<float>("pnorm");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["verbose"] = GetAttr<int>("verbose");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["nlopt"] = GetAttr<int>("nlopt");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["lsa"] = GetAttr<float>("lsa");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["lsb"] = GetAttr<float>("lsb");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["alpha"] = GetAttr<float>("ftalpha");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["m"] = GetAttr<size_t>("ftm");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["ftiter"] = GetAttr<int>("ftiter");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["fteps"] = GetAttr<float>("fteps");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["parallel_linesearch"] = GetAttr<bool>("parallel_linesearch");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["ft"] = GetAttr<int>("ft");
	} catch (const TinyXMLQueryException&) {}
	try {
		_tf = GetAttr<size_t>("frame_duration");
	} catch (const TinyXMLQueryException&) {}
	try {
		ft_params["tv1"] = RHSList<size_t>("tv");
	} catch (const TinyXMLQueryException&) {}
	try {
		_ta = GetAttr<float>("ta");
	} catch (const TinyXMLQueryException&) {
		_ta = wspace.PGet<float>("TA");
	}



	printf ("... done.\n\n");
	m_initialised = true;
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
	Vector<size_t> order(4); order[0]=0; order[1]=2; order[2]=3; order[3]=1;

    size_t nline, nv, nt, nx, nz, nc;
    Vector<size_t> n = size(data);
    nx = n[0]; nv = n[1]; nz = n[2]; nc = n[3];

    // Timing
	nt = ceil(_ta/_tf);
	nv = floor (nv/nt);

	std::cout << "  Data acquired within " << _ta << "s" << std::endl;
	std::cout << "  Sorting data in " << nt << " contrast gates with "
			  << nv << " views each." << std::endl;

	// Clip over the time data, kspace and signal
	data = data(CR(),CR(0,nt*nv-1),CR(),CR());
	kspace = kspace(CR(),CR(),CR(0,nt*nv-1));

	std::cout << "  Incoming    " << std::endl;
	std::cout << "    data:     " << size(data) << std::endl;
	std::cout << "    k-space:  " << size(kspace) << std::endl;

	// Contrasts sorting
	data   = resize(data,nx,nv,nt,nz,nc);
	kspace = resize(kspace,size(kspace,0),nx,nv,nt);

	std::cout << "  Reshaped:   " << std::endl;
	std::cout << "    data:     " << size(data) << std::endl;
	std::cout << "    k-space:  " << size(kspace) << std::endl;

	data   = resize(data,nx*nv,nt,nz,nc);
	data   = permute (data,order);
	kspace = resize(kspace,size(kspace,0),nx*nv,nt);
	weights = zeros<float>(nx*nv,1);
    weights (R( 0,nx/2-1),0) = linspace<float>(1.,0.01,nx/2);
    weights (R(nx/2,nx-1),0) = linspace<float>(0.01,1.,nx/2);
    for (auto i = weights.Begin()+nx; i < weights.End(); i += nx)
        std::copy (weights.Begin(), weights.Begin()+nx, i);

	std::cout << "  Reshaped and permuted:    " << std::endl;
	std::cout << "    data:     " << size(data) << std::endl;
	std::cout << "    k-space:  " << size(kspace) << std::endl;
	std::cout << "    weights:  " << size(weights) << std::endl;

	// FT operator
	Matrix<cxfl>& sensitivities = Get<cxfl>("sensitivities");
	ft_params["sensitivities"] = sensitivities;
	Vector<size_t> image_size = size(sensitivities);
	image_size.pop_back(); image_size.pop_back();
	ft_params["imsz"] = image_size;
	ft_params["dim4"] = (int) nt;
	ft_params["nk"]   = size(data,0);

	Matrix<cxfl> im_xd (image_size[0], image_size[1], nz, nt);
	for (size_t i = 0; i < 4; ++i) {
        ft_params["sensitivities"] = squeeze(sensitivities(CR(),CR(),CR(i),CR()));
        CS_XSENSE<cxfl> ft(ft_params);
        ft.KSpace(kspace);
        ft.Weights (weights);
		im_xd(R(),R(),R(i),R()) = ft ->* squeeze(data(CR(),CR(i),CR(),CR()));
	}

    Add ("im_xd", im_xd);

    return codeare::OK;

}

GRASP::~GRASP() {}


codeare::error_code GRASP::Finalise() {return codeare::OK;}


// the class facxflories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new GRASP;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}


