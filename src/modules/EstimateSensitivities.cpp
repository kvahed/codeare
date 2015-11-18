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

#include "EstimateSensitivities.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Access.hpp"
#include "IOContext.hpp"
#include "Interpolate.hpp"
#include "Math.hpp"

#ifndef NYUGA
	#define NYUGA 111.246117975
#endif

using namespace RRStrategy;

#define clear(X) Free("X")

EstimateSensitivities::EstimateSensitivities () : m_dft_3rd_dim(false), m_test_case(false), m_dim(2) {}


EstimateSensitivities::~EstimateSensitivities () {
	this->Finalise();
}


codeare::error_code EstimateSensitivities::Finalise () {
	return codeare::OK;
}



codeare::error_code EstimateSensitivities::Init () {

	codeare::error_code error = codeare::OK; 
	m_initialised             = false;

	Params p;
	p["alpha"]  = GetAttr<float>("ftalpha");
	p["m"]      = GetAttr<size_t>("ftm");
	p["nk"]     = GetAttr<size_t>("nk")*GetAttr<size_t>("nl");
	p["3rd_dim_cart"] = GetAttr<bool>("cart_3rd_dim");
	p["imsz"]   = GetList<size_t>("image_size");
	p["maxit"]   = GetList<size_t>("ftiter");

	ft = NFFT<cxfl> (p);

	m_initialised = true;

	Matrix<cxfl> sensitivities;
	Add ("sensitivities", sensitivities);

	return error;

}

codeare::error_code EstimateSensitivities::Prepare () {

	codeare::error_code error = codeare::OK;

    Matrix<float>& kspace = Get<float>("sync");
    kspace = squeeze(kspace(CR(),CR(0)));
    kspace = permute(resize(kspace,2,numel(kspace)/2),1,0);

	std::cout << ft << std::endl;

	return error;

}


codeare::error_code EstimateSensitivities::Process () {

    Matrix<float> sos_all, xi, XI;
	Matrix<cxfl>& data = Get<cxfl>("meas");
    Matrix<float>& kspace = Get<float>("sync");
    size_t nk, nl, nc, nv;

    // Normalise data according to sampling time
    // Get rid of singleton dimensions
    std::cout << "    * Get rid of singleton dimensions" << std::endl;
    data *= 1e6;
    data.Squeeze();
    nk = size(data,0);
    nc = size(data,1);
    nl = size(data,2);
    nv = size(data,3);   

    // Forming outgoing data
    Vector<size_t> image_size = GetList<size_t>("image_size");
    image_size.push_back(nc);
	Matrix<cxfl> sensitivities (image_size);
    std::cout << "    * Reconstructing sensitivities to " << size(sensitivities) << std::endl;

    // Interpolate k-space and normalize
    std::cout << "    * Interpolate k-space from gradient raster to sampling int" << std::endl;
    xi = linspace<float>(0.0,1.0,size(kspace,0));
    XI = linspace<float>(0.0,1.0,size(kspace,0)*10);
    kspace  = interp1(xi, kspace, XI);
    kspace  = kspace(CR(0,size(data,0)-1),CR());

    // Create golden angle k-space
    std::cout << "    * Calculate golden angles" << std::endl;
    Matrix<cxfl> ktmp = repmat(ccomplex(kspace(CR(),CR(0)),kspace(CR(),CR(1))),1,nv);
    Matrix<cxfl> ga = repmat(resize(cpolar(1.0f,exp(linspace<float>(0,nv-1,nv)*NYUGA)),1,nv),nk,1);
    ktmp *= ga;
    ktmp /= 2.*mmax(abs(ktmp));
    ktmp  = resize(ktmp,nk*nv,1);
    kspace = Matrix<float>(nk*nv,2);
    kspace (R(),0) = real(ktmp);
    kspace (R(),1) = imag(ktmp);
    kspace = permute (kspace,1,0);

    // Prepare FT operator
    std::cout << "    * Prepare FT operator" << std::endl;
    ft.KSpace(kspace);
    ft.Weights(ones<float>(nk*nv,1));

    // Permute dimensions from col,ch,lin,vol to col,vol,lin,ch
    std::cout << "    * Permute dimensions from col,ch,lin,vol to col,vol,lin,ch" << std::endl;
    
    Matrix<cxfl> tmp (nk, nv, nl, nc);
    for (size_t vol = 0; vol < nv; ++vol)
        for (size_t lin = 0; lin < nl; ++lin)
            for (size_t cha = 0; cha < nc; ++cha)
                std::copy (&data(0,cha,lin,vol),&data(0,cha,lin,vol)+nk,&tmp(0,vol,lin,cha));
    data = Matrix<cxfl>();

    std::cout << "    * Reconstruct coils" << std::endl;
    // Reconstruct all coils
    for (size_t i = 0; i < 1; ++i) {
        sensitivities(R(),R(),R(),R(i)) = ft ->* tmp(CR(),CR(),CR(),CR(i));
        std::cout << "      * 1" << std::endl;
    }

    // Build sum of squares of all coils
    std::cout << "    * Build SOS of sensitivities" << std::endl;
    sos_all = sos(sensitivities);
    
    // Devide images by sos
    std::cout << "    * Devide by SOS" << std::endl;
/*    for (size_t i = 0; i < size(data,4); ++i)
      sensitivities(R(),R(),R(),R(i)) /= sos_all;*/

    // Normalise
    std::cout << "    * Normalise sensitivities" << std::endl;
//    sensitivities /=max(abs(sensitivities));

    Add ("sensitivities", sensitivities);
//    Add ("sos_all", sos_all);
    Add ("kspace", kspace);
	return codeare::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new EstimateSensitivities;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

