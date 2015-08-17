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
#ifdef HAVE_NFFT3
	#include "NCSENSE.hpp"
#endif
#include "Algos.hpp"
#include "Creators.hpp"


using namespace RRStrategy;

Matrix<cxfl> ffdbx, ffdbg, ttdbx, ttdbg, wx, wdx;

template<class T> inline static typename TypeTraits<T>::RT Obj (
    const Matrix<T>& ffdbx, const Matrix<T>& ffdbg, const Matrix<T>& data,
    const typename TypeTraits<T>::RT t) {

	Matrix<T> om = ffdbx;
	if (t > 0.0)
		om += t * ffdbg;
	om -= data;
	return real(om.dotc(om));

}


template<class T> inline static typename TypeTraits<T>::RT TV (
    const Matrix<T>& ttdbx, const Matrix<T>& ttdbg, const typename TypeTraits<T>::RT t,
    const CSParam& cgp) {

	typename TypeTraits<T>::RT o = 0.0, p = 0.5*cgp.pnorm;
	Matrix<T> om = ttdbx;
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
template<class T> inline static typename TypeTraits<T>::RT XFM (
    const Matrix<T>& x, const Matrix<T>& g, const typename TypeTraits<T>::RT t,
    const CSParam& cgp) {
    
	typename TypeTraits<T>::RT o = 0.0, p = 0.5*cgp.pnorm;
	Matrix<T> om = x;
	if (t > 0.0)
		om += t * g;
	om *= conj(om);
	om += cgp.l1;
	om ^= p;
	for (size_t i = 0; i < om.Size(); i++)
		o += om[i].real();
	return cgp.xfmw * o;

}


template<class T> inline typename TypeTraits<T>::RT CompressedSensing::f (
    const Matrix<T>& x, const Matrix<T>& dx, const typename TypeTraits<T>::RT& t,
    typename TypeTraits<T>::RT& rmse) {
    
	typename TypeTraits<T>::RT obj = Obj (ffdbx, ffdbg, data, t);
	rmse = sqrt(obj/m_ndnz);
	if (m_csparam.tvw)
		obj += TV (ttdbx, ttdbg, t, m_csparam);
	if (m_csparam.xfmw)
		obj += XFM (x, dx, t, m_csparam);
	return obj;

}


/**
 * @brief Compute gradient of the data consistency
 */
template<class T> inline static Matrix<T> dObj (
    const Matrix<T>& x, const Matrix<T>& wx, const Matrix<T>& data, const CSParam& cgp) {
    
	FT<T>& ft = *cgp.ft;
    DWT<T>& dwt = *cgp.dwt;
	return (2.0 * (dwt * (ft ->* ((ft * wx) - data))));

}


/**
 * @brief Compute gradient of L1-transform operator
 *
 * @param  x   X
 * @param  cgp CG parameters
 * @return     The gradient
 */
template <class T> inline static Matrix<T> dXFM (
    const Matrix<T>& x, const CSParam& cgp) {
    
	typename TypeTraits<T>::RT pn = 0.5*cgp.pnorm-1.0, l1 = cgp.l1, xfm = cgp.xfmw;
	return xfm * (x * ((x * conj(x) + l1) ^ pn));

}


/**
 * @brief Compute gradient of the total variation operator
 *
 * @param  x   Image space original
 * @param  wx  Image space perturbance
 * @param  cgp Parameters
 */
template<class T> inline static Matrix<T> dTV (
    const Matrix<T>& x, const Matrix<T>& wx, const CSParam& cgp) {

	DWT<T>& dwt = *cgp.dwt;
	TVOP<T>& tvt = *cgp.tvt;
	typename TypeTraits<T>::RT p = 0.5*cgp.pnorm-1.0;
	Matrix<T> dx, g;
	dx = tvt * wx;
	g  = dx * conj(dx);
	g += cgp.l1;
	g ^= p;
	g *= dx;
	g *= cgp.pnorm;
	g  = dwt * (tvt->*g);
	return (cgp.tvw * g);

}


template<class T> inline Matrix<T> CompressedSensing::df (const Matrix<T>& x) {

	DWT<T>& dwt = *m_csparam.dwt;
    wx = dwt->*x;
	Matrix<T> g = dObj (x, wx, data, m_csparam);
	if (m_csparam.xfmw)
		g += dXFM (x, m_csparam);
	if (m_csparam.tvw)
		g += dTV  (x, wx, m_csparam);
	return g;

}

template<class T> inline void CompressedSensing::Update (const Matrix<T>& dx) {
	DWT<cxfl>&  dwt = *m_csparam.dwt;
	FT<cxfl>&   ft  = *m_csparam.ft;
	TVOP<cxfl>& tvt = *m_csparam.tvt;
    wdx =  dwt->*dx;
    ffdbx = ft * wx;
    ffdbg = ft * wdx;
    if (m_csparam.tvw) {
        ttdbx = tvt * wx;
        ttdbg = tvt * wdx;
    }
}

template<class T> inline void CompressedSensing::NLCG (Matrix<T>& x) {
    
	typename TypeTraits<T>::RT t0  = 1.0, t = 1.0, z = 0., xn = norm(x), rmse, bk, f0, f1, dxn;
	Matrix<T> dx, g0, g1;
    
	g0  = df (x);
	dx  = -g0;
    
	for (size_t k = 0; k < (size_t)m_csparam.cgiter; k++) {
        
        Update(dx);
        
		t = t0;
        
		f0 = f (x, dx, z, rmse);
        
		int i = 0;
		while (i < m_csparam.lsiter) {
			t *= m_csparam.lsb;
			f1 = f (x, dx, t, rmse);
			if (f1 <= f0 - (m_csparam.lsa * t * abs(g0.dotc(dx))))
				break;
			++i;
		}
        
		printf (ofstr.c_str(), k, rmse, i); fflush (stdout);
        
		if (i == m_csparam.lsiter) {
			printf ("Reached max line search, exiting... \n");
			return;
		}
        
		if      (i > 2) t0 *= m_csparam.lsb;
		else if (i < 1) t0 /= m_csparam.lsb;
        
		// Update image
		x  += dx * t;
        
		// CG computation
		g1  =  df (x);
		bk  =  real(g1.dotc(g1)) / real(g0.dotc(g0));
		g0  =  g1;
		dx  = -g1 + dx * bk;
		dxn =  norm(dx)/xn;
        
		printf ("dxnrm: %0.4f\n", dxn);
		if (dxn < m_csparam.cgconv)
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
#ifdef HAVE_NFFT3
			ft_params["epsilon"] = RHSAttribute<double>("fteps");
			ft_params["alpha"]   = RHSAttribute<double>("ftalpha");
			ft_params["maxit"]   = RHSAttribute<size_t>("ftiter");
			ft_params["m"]       = RHSAttribute<size_t>("ftm");
	        ft_params["nk"]      = RHSAttribute<size_t>("ftnk");
	        ft_params["3rd_dim_cart"] = RHSAttribute<bool>("cart_3rd_dim");
		    ft_params["threads"] = RHSAttribute<int>("threads");
	        ft_params["imsz"]    = m_image_size;
	        m_csparam.ft = (FT<cxfl>*) new NFFT<cxfl> (ft_params);
#else
			printf("**ERROR - CompressedSensing: NUFFT support not available.");
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

	data = (m_test_case) ?
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
    m_ndnz = (float)nnz(data);
	im_dc  = dwt * im_dc;

	printf ("  Running %i NLCG iterations ... \n", m_csiter); fflush(stdout);

	for (size_t i = 0; i < (size_t)m_csiter; i++) {
		NLCG (im_dc);
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


