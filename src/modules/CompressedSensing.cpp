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

#include "CompressedSensing.hpp"
#include "Toolbox.hpp"
#include "DFT.hpp"
#ifdef HAVE_NFFT
	#include "NCSENSE.hpp"
#endif
#include "Algos.hpp"
#include "Creators.hpp"

using namespace RRStrategy;

codeare::error_code CompressedSensing::Init () {

	printf ("Intialising CompressedSensing ...\n");

	Params ft_params;

	for (size_t i = 0; i < 3; i++)
		m_N[i] = 1;

	Attribute ("tvw",     &m_csparam.tvw);
	Attribute ("xfmw",    &m_csparam.xfmw);
	Attribute ("l1",      &m_csparam.l1);
	Attribute ("pnorm",   &m_csparam.pnorm);
	Attribute ("verbose", &m_verbose);
    m_image_size   = RHSList<size_t>("ftdims");

	printf ("  Geometry: " JL_SIZE_T_SPECIFIER "D (" JL_SIZE_T_SPECIFIER ","
            JL_SIZE_T_SPECIFIER "," JL_SIZE_T_SPECIFIER ")\n", m_image_size.size(),
            m_image_size[0], m_image_size[1], (m_image_size.size()==3) ? m_image_size[2] : 1);

	Attribute ("test_case", &m_test_case);
	Attribute ("noise",     &m_noise);

    printf ("  Weights: TV(%.2e) XF(%.2e) L1(%.2e)\n", m_csparam.tvw, m_csparam.xfmw, m_csparam.l1);
    printf ("  Pnorm: %.2e\n", m_csparam.pnorm);
	
    Attribute ("ft", &m_ft_type);

	printf ("  FFT class: ");
	switch (m_ft_type)
		{
		case 0:
			printf ("%s", "DFT");
			ft_params["dims"] = m_image_size;
			ft_params["mask"] = Get<float>("mask");
		    ft_params["threads"] = RHSAttribute<int>("threads");
			m_csparam.ft = (FT<cxfl>*) new DFT<cxfl> (ft_params);
			break;
		case 1:
			printf ("%s", "SENSE");
			assert(false);
			break;
		case 2:
			printf ("%s", "NUFFT");
#ifdef HAVE_NFFT
			ft_params["epsilon"] = RHSAttribute<double>("fteps");
			ft_params["alpha"]   = RHSAttribute<double>("ftalpha");
			ft_params["maxit"]   = RHSAttribute<size_t>("ftiter");
			ft_params["m"]       = RHSAttribute<size_t>("ftm");
	        ft_params["nk"]      = RHSAttribute<size_t>("ftnk");
	        ft_params["imsz"]    = m_image_size;
	        m_csparam.ft = (FT<cxfl>*) new NFFT<cxfl> (ft_params);
#else
			printf("**ERROR - CompressedSensing: NUFFT support not available.");
			assert(false);
#endif
			break;
		case 3:
			printf ("%s", "NCSENSE");
#ifdef HAVE_NFFT
			ft_params["sensitivities"] = Get<cxfl>("sensitivities");
            ft_params["nk"]           = (size_t) RHSAttribute<int>("nk");
			ft_params["weights_name"] = std::string("weights");
		    ft_params["ftiter"]       = (size_t) RHSAttribute<int>("ftmaxit");
		    ft_params["fteps"]        = RHSAttribute<double>("fteps");
		    ft_params["cgiter"]       = (size_t) RHSAttribute<int>("cgmaxit");
		    ft_params["cgeps"]        = RHSAttribute<double>("cgeps");
		    ft_params["lambda"]       = RHSAttribute<double>("lambda");
		    ft_params["threads"]      = RHSAttribute<int>("threads");
			m_csparam.ft = (FT<cxfl>*) new NCSENSE<cxfl> (ft_params);
#else
			printf("**ERROR - CompressedSensing: NUFFT support not available.");
			assert(false);
#endif
			break;
		default:
			printf ("No FT strategy defined");
			assert (false);
			break;
		}
	printf ("\n");

	Attribute ("csiter",    &m_csiter);
	Attribute ("wl_family", &m_wf);
	Attribute ("wl_member", &m_wm);
	printf ("  DWT(%i,%i)", m_wf, m_wm);
	
	if (m_wf < -1 || m_wf > 5)
		m_wf = -1;

	Attribute ("cgconv", &m_csparam.cgconv);
	Attribute ("cgiter", &m_csparam.cgiter);
	Attribute ("lsiter", &m_csparam.lsiter);
	Attribute ("lsa",    &m_csparam.lsa);
	Attribute ("lsb",    &m_csparam.lsb);
    printf ("  Iterations: CS(%i) CG(%i) LS(%i)\n", m_csiter, m_csparam.cgiter, m_csparam.lsiter);
	printf ("  Conv: CG(%.4f)\n", m_csparam.cgconv);
	printf ("  LS brackets: lsa(%.2e) lsb(%.2e)", m_csparam.lsa, m_csparam.lsb);

	m_initialised = true;
	printf ("... done.\n\n");

	return codeare::OK;

}


codeare::error_code CompressedSensing::Prepare () {

	codeare::error_code error = codeare::OK;

	FT<cxfl>& dft = *m_csparam.ft;

	if (m_ft_type == 2 || m_ft_type == 3) {
		dft.KSpace (Get<float>("kspace"));
		dft.Weights (Get<float>("weights"));
	} else {
		dft.Mask (Get<float>("mask"));
	}

	Free ("weights");
	Free ("kspace");

	m_initialised = true;

	return error;

}

codeare::error_code CompressedSensing::Process () {

	float ma;

	FT<cxfl>& dft = *m_csparam.ft;

	Matrix<cxfl> data  = m_test_case ?
		dft * phantom<cxfl>(m_image_size[0]) : Get<cxfl>("data");

	if (m_noise > 0.)
		data += m_noise * randn<cxfl>(size(data));

	Matrix<float>& pdf   = Get<float>  ("pdf" );
	Matrix<cxfl>&  pc    = Get<cxfl>   ("pc");
    Matrix<cxfl> im_dc;

    m_csparam.dwt = new DWT <cxfl> (m_image_size[0], (wlfamily) m_wf, m_wm);
	m_csparam.tvt = new TVOP ();

	DWT<cxfl>& dwt = *m_csparam.dwt;
    std::vector< Matrix<cxfl> > vc;
	
	im_dc  = data;
	if (m_ft_type != 2 && m_ft_type != 3)
		im_dc /= pdf;
	im_dc  = dft ->* im_dc;

	ma       = m_max(abs(im_dc));

	if (m_verbose)
		vc.push_back(im_dc);

	im_dc /= ma;
	data  /= ma;
	im_dc  = dwt * im_dc;

	printf ("  Running %i NLCG iterations ... \n", m_csiter); fflush(stdout);

	for (size_t i = 0; i < (size_t)m_csiter; i++) {
		NLCG (im_dc, data, m_csparam);
		if (m_verbose)
			vc.push_back(dwt ->* im_dc*ma);
	}

    if (m_verbose) {
        size_t cpsz = numel(im_dc);
        im_dc = zeros<cxfl> (size(im_dc,0), size(im_dc,1), (m_dim == 3) ?
        		size(im_dc,2) : 1, vc.size());
        for (size_t i = 0; i < vc.size(); i++)
            memcpy (&im_dc[i*cpsz], &(vc[i][0]), cpsz*sizeof(cxfl));

    } else
        im_dc = dwt ->* im_dc * ma;

    Add ("im_dc", im_dc);

    return codeare::OK;

}


CompressedSensing::CompressedSensing() :
	m_wm(0), m_csiter(0), m_wf(0), m_dim(0), m_verbose(0), m_ft_type(0),
	m_noise(0.), m_test_case(0) {}


CompressedSensing::~CompressedSensing() {}


codeare::error_code
CompressedSensing::Finalise() {return codeare::OK;}


// the class facxflories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CompressedSensing;
}


extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}


