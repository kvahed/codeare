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

	Params ft_params;

	for (size_t i = 0; i < 3; i++)
		m_N[i] = 1;

	Attribute ("verbose", &m_verbose);
    m_image_size   = RHSList<size_t>("ftdims");
    m_dim = m_image_size.size();

	printf ("  Geometry: " JL_SIZE_T_SPECIFIER "D (" JL_SIZE_T_SPECIFIER ","
            JL_SIZE_T_SPECIFIER "," JL_SIZE_T_SPECIFIER ")\n", (size_t)m_dim,
            m_image_size[0], m_image_size[1], (m_image_size.size()==3) ? m_image_size[2] : 1);

	Attribute ("test_case", &m_test_case);
	Attribute ("noise",     &m_noise);
	Attribute ("dim4",     &m_dim4);
	Attribute ("dim5",     &m_dim5);

/*    printf ("  Weights: TV(%.2e) XF(%.2e) L1(%.2e)\n", m_csparam.tvw, m_csparam.xfmw, m_csparam.l1);
      printf ("  Pnorm: %.2e\n", m_csparam.pnorm);*/
	
    Attribute ("ft", &m_ft_type);

	printf ("  FFT class: ");
	switch (m_ft_type)
		{
		case 0:
			printf ("%s", "ft");
			ft_params["dims"] = m_image_size;
			ft_params["mask"] = Get<float>("mask");
		    ft_params["threads"] = RHSAttribute<int>("threads");
			break;
		case 1:
			printf ("%s", "SENSE");
			assert(false);
			break;
		case 2:
			printf ("%s", "NUFFT");
#ifdef HAVE_NFFT3
			ft_params["epsilon"] = RHSAttribute<float>("fteps");
			ft_params["alpha"]   = RHSAttribute<float>("ftalpha");
			ft_params["maxit"]   = RHSAttribute<size_t>("ftiter");
			ft_params["m"]       = RHSAttribute<size_t>("ftm");
	        ft_params["nk"]      = RHSAttribute<size_t>("ftnk");
	        ft_params["3rd_dim_cart"] = RHSAttribute<bool>("cart_3rd_dim");
		    ft_params["threads"] = RHSAttribute<int>("threads");
	        ft_params["imsz"]    = m_image_size;
#else
			printf("**ERROR - XDGRASP: NUFFT support not available.");
			assert(false);
#endif
			break;
		case 3:
			printf ("%s", "NCSENSE");
#ifdef HAVE_NFFT3
			ft_params["sensitivities"] = Get<cxfl>("sensitivities");
            ft_params["nk"]           = (size_t) RHSAttribute<int>("nk");
			ft_params["weights_name"] = std::string("weights");
		    ft_params["ftiter"]       = (size_t) RHSAttribute<int>("ftmaxit");
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
#else
			printf("**ERROR - XDGRASP: NUFFT support not available.");
			assert(false);
#endif
			break;
		default:
			printf ("No FT strategy defined");
			assert (false);
			break;
		}
	printf ("\n");

    ft_params["imsz"]  = m_image_size;
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

    //Aaaargh!
    ft_params["wl_family"] = RHSAttribute<int>("wl_family");    
    ft_params["wl_member"] = RHSAttribute<int>("wl_member");
    ft_params["csiter"] = RHSAttribute<int>("csiter");
    ft_params["nliter"] = RHSAttribute<int>("cgiter");
    ft_params["cgconv"] = RHSAttribute<float>("cgconv");
    ft_params["lsiter"] = RHSAttribute<int>("lsiter");
    ft_params["ft"] = RHSAttribute<int>("ft");
    ft_params["dim4"] = m_dim4;
    ft_params["dim5"] = m_dim5;

    csx = new CS_XSENSE<cxfl>(ft_params);


	m_initialised = true;
	std::cout << *csx << std::endl;

	printf ("... done.\n\n");

	return codeare::OK;

}


codeare::error_code XDGRASP::Prepare () {

	codeare::error_code error = codeare::OK;

	FT<cxfl>& ft = *csx;

	if (m_ft_type == 2 || m_ft_type == 3) {
        ft.KSpace(Get<float>("kspace"));
		ft.Weights (Get<float>("weights"));
        Free ("weights");
        Free ("kspace");
    } else {
		ft.Mask (Get<float>("mask"));
    }

	m_initialised = true;

	return error;

}

codeare::error_code XDGRASP::Process () {

    m_image_size.push_back (m_dim4);
    m_image_size.push_back (m_dim5);
    
    Matrix<cxfl> im_xd (m_image_size);
    Matrix<cxfl>& data = Get<cxfl>("data");
    FT<cxfl>& ft = *csx;

    im_xd = ft ->* data;

    Add ("im_xd", im_xd);

    return codeare::OK;

}


XDGRASP::XDGRASP() :
	m_wm(0), m_csiter(0), m_wf(0), m_dim(0), m_verbose(0), m_ft_type(0),
	m_noise(0.), m_test_case(0) {}


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


