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
#ifdef HAVE_NFFT
	#include "NCSENSE.hpp"
#endif
#include "Algos.hpp"
#include "Creators.hpp"


using namespace RRStrategy;

inline static float
Obj (const Matrix<cxfl>& ffdbx, const Matrix<cxfl>& ffdbg,
	 const Matrix<cxfl>&  data, const float             t) {

	Matrix<cxfl> om = ffdbx;

	if (t > 0.0)
		om += t * ffdbg;
	om -= data;

	return real(om.dotc(om));

}


inline static float
ObjTV (const Matrix<cxfl>& ttdbx, const Matrix<cxfl>& ttdbg,
	   const float             t, const CSParam&        cgp) {

	float o = 0.0, p = 0.5*cgp.pnorm;
	Matrix<cxfl> om = ttdbx;

	if (t > 0.0)
		om += t * ttdbg;
	om *= conj(om);
	om += cgp.l1;
	om ^= p;

	for (size_t i = 0; i < om.Size(); i++)
		o += real(om[i]);

	return cgp.tvw * o;

}


/**
 *
 */
inline static float
ObjXFM (const Matrix<cxfl>& x, const Matrix<cxfl>& g,
		const float         t, const CSParam&    cgp) {

	float o = 0.0, p = 0.5*cgp.pnorm;
	Matrix<cxfl> om = x;

	if (t > 0.0)
		om += t * g;
	om *= conj(om);
	om += cgp.l1;
	om ^= p;

	for (size_t i = 0; i < om.Size(); i++)
		o += om[i].real();

	return cgp.xfmw * o;

}


static float
Objective (const Matrix<cxfl>& ffdbx, const Matrix<cxfl>& ffdbg,
		   const Matrix<cxfl>& ttdbx, const Matrix<cxfl>& ttdbg,
		   const Matrix<cxfl>&     x, const Matrix<cxfl>&     g,
		   const Matrix<cxfl>&  data, const float             t,
				 float&         rmse, const CSParam&        cgp) {

	float obj = Obj (ffdbx, ffdbg, data, t);

	rmse = sqrt(obj/(float)nnz(data));

	if (cgp.tvw)
		obj += ObjTV (ttdbx, ttdbg, t, cgp);
	if (cgp.xfmw)
		obj += ObjXFM (x, g, t, cgp);

	return obj;

}


/**
 * @brief Compute gradient of the data consistency
 */
static Matrix<cxfl>
GradObj (const Matrix<cxfl>& x, const Matrix<cxfl>& wx,
		 const Matrix<cxfl>& data, const CSParam& cgp) {

	FT<cxfl>& ft = *(cgp.ft);
	DWT<cxfl>& dwt = *(cgp.dwt);

	return (2.0 * (dwt * (ft ->* ((ft * wx) - data))));

}


/**
 * @brief Compute gradient of L1-transform operator
 *
 * @param  x   X
 * @param  cgp CG parameters
 * @return     The gradient
 */
template <class T> inline static Matrix<T>
GradXFM   (const Matrix<T>& x, const CSParam& cgp) {

	float pn = 0.5*cgp.pnorm-1.0, l1 = cgp.l1, xfm = cgp.xfmw;
	return xfm * (x * ((x * conj(x) + l1) ^ pn));

}


/**
 * @brief Compute gradient of the total variation operator
 *
 * @param  x   Image space original
 * @param  wx  Image space perturbance
 * @param  cgp Parameters
 */
Matrix<cxfl>
GradTV    (const Matrix<cxfl>& x, const Matrix<cxfl>& wx, const CSParam& cgp) {

	DWT<cxfl>& dwt = *cgp.dwt;
	TVOP<cxfl>& tvt = *cgp.tvt;
	float p   = 0.5*cgp.pnorm-1.0;
	Matrix<cxfl> dx, g;


	dx = tvt * wx;
	g  = dx * conj(dx);
	g += cgp.l1;
	g ^= p;
	g *= dx;
	g *= cgp.pnorm;
	g  = dwt * (tvt->*g);

	return (cgp.tvw * g);

}


Matrix<cxfl> Gradient (const Matrix<cxfl>& x, const Matrix<cxfl>& wx,
		const Matrix<cxfl>& data, const CSParam& cgp) {

	Matrix<cxfl> g = GradObj (x, wx, data, cgp);

	if (cgp.xfmw)
		g += GradXFM (x, cgp);
	if (cgp.tvw)
		g += GradTV  (x, wx, cgp);
	return g;

}


void NLCG (Matrix<cxfl>& x, const Matrix<cxfl>& data, const CSParam& cgp) {


	float     t0  = 1.0, t = 1.0, z = 0., xn = norm(x), rmse, bk, f0, f1, dxn;

	Matrix<cxfl> g0, g1, dx, ffdbx, ffdbg, ttdbx, ttdbg, wx, wdx;

	DWT<cxfl>& dwt = *cgp.dwt;
	FT<cxfl>&  ft  = *cgp.ft;
	TVOP<cxfl>&      tvt = *cgp.tvt;

	wx  = dwt->*x;

	g0 = Gradient (x, wx, data, cgp);
	dx = -g0;
	wdx = dwt->*dx;

	for (size_t k = 0; k < (size_t)cgp.cgiter; k++) {

		t = t0;

		ffdbx = ft * wx;
		ffdbg = ft * wdx;

		if (cgp.tvw) {
			ttdbx = tvt * wx;
			ttdbg = tvt * wdx;
		}

		f0 = Objective (ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data, z, rmse, cgp);

		int i = 0;
		while (i < cgp.lsiter) {

			t *= cgp.lsb;
			f1 = Objective(ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data, t, rmse, cgp);
			if (f1 <= f0 - (cgp.lsa * t * abs(g0.dotc(dx))))
				break;
			++i;
		}

		printf (ofstr.c_str(), k, rmse, i); fflush (stdout);

		if (i == cgp.lsiter) {
			printf ("Reached max line search, exiting... \n");
			return;
		}

		if      (i > 2) t0 *= cgp.lsb;
		else if (i < 1) t0 /= cgp.lsb;

		// Update image
		x  += (dx * t);
		wx  = dwt->*x;

		// CG computation
		g1  =  Gradient (x, wx, data, cgp);
		bk  =  real(g1.dotc(g1)) / real(g0.dotc(g0));
		g0  =  g1;
		dx  = -g1 + dx * bk;
		wdx =  dwt->*dx;
		dxn =  norm(dx)/xn;

		printf ("dxnrm: %0.4f\n", dxn);
		if (dxn < cgp.cgconv)
			break;

	}


}


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
    m_dim = m_image_size.size();

	printf ("  Geometry: " JL_SIZE_T_SPECIFIER "D (" JL_SIZE_T_SPECIFIER ","
            JL_SIZE_T_SPECIFIER "," JL_SIZE_T_SPECIFIER ")\n", (size_t)m_dim,
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
			printf ("%s", "ft");
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
	        ft_params["3rd_dim_cart"] = RHSAttribute<bool>("cart_3rd_dim");
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
	        ft_params["3rd_dim_cart"] = RHSAttribute<bool>("cart_3rd_dim");
		    ft_params["verbose"]      = 0;
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

	FT<cxfl>& ft = *m_csparam.ft;

	if (m_ft_type == 2 || m_ft_type == 3) {
		ft.KSpace (Get<float>("kspace"));
		ft.Weights (Get<float>("weights"));
	} else {
		ft.Mask (Get<float>("mask"));
	}

	std::cout << ft << std::endl;

	Free ("weights");
	Free ("kspace");

	m_initialised = true;

	return error;

}

codeare::error_code CompressedSensing::Process () {

	float ma;
	FT<cxfl>& ft = *m_csparam.ft;

	Matrix<cxfl> data = (m_test_case) ?
		ft * phantom<cxfl>(m_image_size[0], m_image_size[0], (m_dim == 3) ?
				m_image_size[0] : 1) : Get<cxfl> ("data");

	if (m_noise > 0.)
		data += m_noise * randn<cxfl>(size(data));

	Matrix<float>& pdf   = Get<float>("pdf");
	Matrix<cxfl> im_dc;

    m_csparam.dwt = new DWT <cxfl> (m_image_size[0], (wlfamily)m_wf, m_wm);
	m_csparam.tvt = new TVOP<cxfl> ();

	DWT<cxfl>& dwt = *m_csparam.dwt;
    std::vector< Matrix<cxfl> > vc;
	
	im_dc  = data;
	if (m_ft_type != 2 && m_ft_type != 3)
		im_dc /= pdf;
	im_dc  = ft ->* im_dc;

	ma     = max(abs(im_dc));

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
        	std::copy (&vc[i][0], &vc[i][0]+cpsz, &im_dc[i*cpsz]);
    } else
        im_dc = dwt ->* im_dc*ma;

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


