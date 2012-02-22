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

#ifndef __DWT_HPP__
#define __DWT_HPP__

enum wlfamily {
	
	ID = -1,
	WL_DAUBECHIES,
	WL_DAUBECHIES_CENTERED,
	WL_HAAR,
	WL_HAAR_CENTERED,
	WL_BSPLINE,
	WL_BSPLINE_CENTERED

};

#include "Matrix.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_wavelet2d.h>



/**
 * @brief 2D Discrete wavelet transform for Matrix template<br/>(Daubechies wavelets)
 */
class DWT {
	

public:


	/**
	 * @brief Construct 2D Wavelet transform with wavelet class and side length
	 *
	 * @param  sl      Side length
	 * @param  wf      Wavelet family (default none, i.e. ID)
	 * @param  wm      Familty member (default 4)
	 */
	DWT (const size_t& sl, const wlfamily& wf = ID, const size_t& = 4);


	virtual 
	~DWT();


	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<cxfl> 
	Trafo        (const Matrix<cxfl>& m) const ;
	

	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<cxfl> 
	Adjoint      (const Matrix<cxfl>& m) const ;
	

private:
	
	DWT();


	/**
	 * @brief   Transform
	 *
	 * @param   m   To transform
	 * @param   bw  Backward: true, Forward: false
	 * @return      Transform
	 */
	Matrix<cxfl> 
	Transform    (const Matrix<cxfl>& m, const bool& bw) const ;

	wlfamily m_wf;

	size_t  m_sz;
	size_t  m_sl;

	double* m_re;
	double* m_im;
	
	gsl_wavelet_workspace* m_work;
	gsl_wavelet           *m_w;
	
};

#endif
