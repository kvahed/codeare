/*
w *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
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

#include "XDGRASP.hpp"
#include "Toolbox.hpp"
#ifdef HAVE_NFFT3
	#include "NCSENSE.hpp"
#endif
#include "CS_XSENSE.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

using namespace RRStrategy;

codeare::error_code XDGRASP::Init () {

	printf ("Intialising XDGRASP ...\n");

    std::cout << *this << std::endl;

	for (size_t i = 0; i < 3; i++)
		m_N[i] = 1;

	Attribute ("verbose", &m_verbose);

	Attribute ("test_case", &m_test_case);
	Attribute ("noise",     &m_noise);


	ft_params["weights_name"] = std::string("weights");
	ft_params["alpha"]   = RHSAttribute<float>("ftalpha");
	ft_params["maxit"]   = RHSAttribute<size_t>("ftiter");
	ft_params["ftiter"]       = (size_t) RHSAttribute<int>("ftmaxit");
	ft_params["m"]       = RHSAttribute<size_t>("ftm");
	ft_params["fteps"]        = RHSAttribute<float>("fteps");
	ft_params["cgiter"]       = (size_t) RHSAttribute<int>("cgmaxit");
	ft_params["cgeps"]        = RHSAttribute<float>("cgeps");
	ft_params["lambda"]       = RHSAttribute<float>("lambda");
	ft_params["threads"]      = RHSAttribute<int>("threads");
	ft_params["3rd_dim_cart"] = RHSAttribute<bool>("cart_3rd_dim");
	try {
		size_t n  = GetAttr<size_t> ("nmany");
		ft_params["nmany"] = n;
	} catch (const TinyXMLQueryException&) {
		ft_params["nmany"] = 1;
	}

	ft_params["verbose"]      = 0;

    ft_params["nlopt"] = RHSAttribute<int>("nlopt");
    ft_params["tvw1"] = RHSAttribute<float>("tvw1");
    ft_params["tv1"] = RHSList<size_t>("tv1");
    ft_params["tvw2"] = RHSAttribute<float>("tvw2");
    ft_params["tv2"] = RHSList<size_t>("tv2");
    ft_params["xfmw"] = RHSAttribute<float>("xfmw");
    ft_params["l1"] = RHSAttribute<float>("l1");
    ft_params["lsa"] = RHSAttribute<float>("lsa");
    ft_params["lsb"] = RHSAttribute<float>("lsb");
    ft_params["pnorm"] = RHSAttribute<float>("pnorm");    
    ft_params["threads"] = RHSAttribute<int>("threads");
    ft_params["verbose"] = RHSAttribute<int>("verbose");
    ft_params["parallel_linesearch"] = RHSAttribute<bool>("parallel_linesearch");
    ft_params["wl_family"] = RHSAttribute<int>("wl_family");    
    ft_params["wl_member"] = RHSAttribute<int>("wl_member");
    ft_params["csiter"] = RHSAttribute<int>("csiter");
    ft_params["nliter"] = RHSAttribute<int>("cgiter");
    ft_params["cgconv"] = RHSAttribute<float>("cgconv");
    ft_params["lsiter"] = RHSAttribute<int>("lsiter");
    ft_params["ft"] = RHSAttribute<int>("ft");

	try {
		_ntres = GetAttr<size_t>("_ntres");
	} catch (const TinyXMLQueryException&) {}

	try {
		_tf = GetAttr<size_t>("frame_duration");
	} catch (const TinyXMLQueryException&) {}

	m_initialised = true;

	printf ("... done.\n\n");

	return codeare::OK;

}

codeare::error_code XDGRASP::Prepare () {
	codeare::error_code error = codeare::OK;
	return error;
}

codeare::error_code XDGRASP::Process () {

    Matrix<cxfl>& data = Get<cxfl>("meas");
    Matrix<float>& kspace = Get<float>("kspace");
    Matrix<float>& weights = Get<float>("weights");
    Matrix<float>& res_signal = Get<float>("res_signal");
    size_t nline, nv, nt, nx, nz, nc;
    float ta = wspace.PGet<float>("TA");
    Vector<size_t> n = size(data);
    nx = n[0]; nv = n[1]; nz = n[2]; nc = n[3];

    // Timing
	nt = ceil(ta/_tf);
	nline = floor(nv/(_ntres*nt));
	std::cout << "  Data acquired within " << ta << std::endl;
	std::cout << "  Sorting data in " << _ntres << " respiratory and " << nt
			  << " contrast gates with "	<< nline << " views each." <<std::endl;

	// Clip over the time data, kspace and signal
	data = data(CR(),CR(0,nt*_ntres*nline-1),CR(),CR());
	res_signal = res_signal(CR(0,nt*_ntres*nline-1));
	kspace = kspace(CR(),CR(),CR(0,nt*_ntres*nline-1));
	nv = size(data,1);

	std::cout << "  Incoming    " << std::endl;
	std::cout << "    data:     " << size(data) << std::endl;
	std::cout << "    kspace:   " << size(kspace) << std::endl;
	std::cout << "    resp sig: " << size(res_signal) << std::endl;

	// Contrasts sorting
	data = resize(data,nx,nv/nt,nt,nz,nc);
	kspace = resize(kspace,size(kspace,0),nx,nv/nt,nt);
	res_signal = resize(res_signal,nv/nt,nt);

	std::cout << "  Reshaped:    " << std::endl;
	std::cout << "    data:     " << size(data) << std::endl;
	std::cout << "    kspace:   " << size(kspace) << std::endl;
	std::cout << "    resp sig: " << size(res_signal) << std::endl;

	// Respiration sorting
	for (size_t i = 0; i < nt; ++i) {
		Matrix<float> res_gate = res_signal(CR(),CR(i));
		Vector<size_t> index = sort(res_gate);
		data(R(),R(),R(i),R(),R()) = data (CR(),CR(index),CR(i),CR(),CR());
		kspace (R(),R(),R(),R(i)) = kspace (CR(),CR(),CR(index),CR(i));
	}
	data = resize(data,nx*nline,nt,_ntres,nz,nc);

	std::cout << "  Reshaped and permuted:    " << std::endl;
	Vector<size_t> order(5); order[0]=0; order[1]=3; order[2]=4; order[3]=1; order[4]=2;
	data = permute (data,order);
	kspace = resize(kspace,size(kspace,0),nx*nline,nt,_ntres);
	weights = zeros<float>(nx*nline,1);
    weights (R( 0,nx/2-1),0) = linspace<float>(1.,1./nx,nx/2);
    weights (R(nx/2,nx-1),0) = linspace<float>(1./nx,1.,nx/2);
    for (auto i = weights.Begin()+nx; i < weights.End(); i += nx)
        std::copy (weights.Begin(), weights.Begin()+nx, i);

	std::cout << "  Reshaped and permuted:    " << std::endl;
	std::cout << "    data:          " << size(data) << std::endl;
	std::cout << "    data:          " << size(data) << std::endl;
	std::cout << "    kspace:        " << size(kspace) << std::endl;
	std::cout << "    weights:       " << size(weights) << std::endl;

	// FT operator
	Matrix<cxfl>& sensitivities = Get<cxfl>("sensitivities");;
	ft_params["sensitivities"] = sensitivities;
	Vector<size_t> image_size = size(sensitivities); image_size.pop_back();
	ft_params["imsz"] = image_size;
	ft_params["dim4"] = (int)nt;
	ft_params["dim5"] = (int)_ntres;
	ft_params["nk"]   = size(data,0);

    CS_XSENSE<cxfl> ft(ft_params);
	std::cout << ft << std::endl;
    ft.KSpace(kspace);
	ft.Weights (weights);

	// 1-3 NuFFT 4,5 TV CS
	Matrix<cxfl> im_xd = ft ->* data;

    Add ("im_xd", im_xd);

    return codeare::OK;

}


XDGRASP::XDGRASP() :
	m_wm(0), m_csiter(0), m_wf(0), m_dim(0), m_verbose(0), m_noise(0.), m_test_case(0),
    _ntres(4), _tf(15.0) {}


XDGRASP::~XDGRASP() {}


codeare::error_code
XDGRASP::Finalise() {return codeare::OK;}


// the class facxflories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new XDGRASP;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}


