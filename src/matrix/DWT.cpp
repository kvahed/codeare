/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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

#include "DWT.hpp"

DWT::DWT (const size_t& sl, const wlfamily& wf, const size_t& wm) {

	// Checks missing !!!

	m_wf = wf;

	if (m_wf > ID) {

		m_sl   = sl;
		m_sz   = sl*sl;

		switch (wf) 
			{
			case WL_DAUBECHIES         : m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies, wm); break;
			case WL_DAUBECHIES_CENTERED: m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies, wm); break;
			case WL_HAAR               : m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies, wm); break;
			case WL_HAAR_CENTERED      : m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies, wm); break;
			case WL_BSPLINE            : m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies, wm); break;
			case WL_BSPLINE_CENTERED   : m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies, wm); break;
			default                    : m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies, wm); break;
			}

		m_work = gsl_wavelet_workspace_alloc (m_sz);
		
		m_re = (double*) malloc (m_sz * sizeof(double));
		m_im = (double*) malloc (m_sz * sizeof(double));

	}

}


DWT::~DWT () {

	if (m_wf != ID) {

		free (m_re);
		free (m_im);
		
		gsl_wavelet_workspace_free (m_work);

	}

}


Matrix<cxfl> 
DWT::Trafo (const Matrix<cxfl>& m) const {

	return Transform (m, false);

}


Matrix<cxfl> 
DWT::Adjoint (const Matrix<cxfl>& m) const {

	return Transform (m, true);

}


Matrix<cxfl> 
DWT::Transform (const Matrix<cxfl>& m, const bool& bw) const {
	
	Matrix<cxfl> res = m;
	
	if (m_wf > ID) {
		
		// Checks missing !!!
		
		for (size_t j = 0; j < m_sl; j++)
			for (size_t i = 0; i < m_sl; i++) {
				m_re [i*m_sl+j] = res.At(j*m_sl+i).real();
				m_im [i*m_sl+j] = res.At(j*m_sl+i).imag();
			}
		
		if (bw) {
			
			if (!(gsl_wavelet2d_nstransform_inverse (m_w, m_re, m_sl, m_sl, m_sl, m_work) == GSL_SUCCESS))
				printf ("Wavelet transform for real part failed\n.");
			if (!(gsl_wavelet2d_nstransform_inverse (m_w, m_im, m_sl, m_sl, m_sl, m_work) == GSL_SUCCESS))
				printf ("Wavelet transform for imaginary part failed\n.");
			
		} else {
			
			if (!(gsl_wavelet2d_nstransform_forward (m_w, m_re, m_sl, m_sl, m_sl, m_work) == GSL_SUCCESS))
				printf ("Wavelet transform for real part failed\n.");
			if (!(gsl_wavelet2d_nstransform_forward (m_w, m_im, m_sl, m_sl, m_sl, m_work) == GSL_SUCCESS))
				printf ("Wavelet transform for imaginary part failed\n.");
		}
		
		for (size_t j = 0; j < m_sl; j++)
			for (size_t i = 0; i < m_sl; i++) 
				res.At(j*m_sl+i) = cxfl(m_re[i*m_sl+j],m_im [i*m_sl+j]);
		
	}
	
	return res;
	
}



template<> Matrix<cxfl> 
DWT::operator* (const Matrix<cxfl>& m) const {

	return Trafo(m);

}


template<> Matrix<cxfl> 
DWT::operator->* (const Matrix<cxfl>& m) const {

	return Adjoint (m);

}

