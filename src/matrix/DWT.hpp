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

/**
 * @brief  Supported wavelet families
 */
enum wlfamily {
	
	ID = -1,                  /**< Identity transform*/
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
	 * @brief    Adjoint transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<cxfl> 
	Adjoint      (const Matrix<cxfl>& m) const ;
	

	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	template <class T> Matrix<T> 
	operator*    (const Matrix<T>& m) const ;
	

	/**
	 * @brief    Adjoint transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	template <class T> Matrix<T> 
	operator->* (const Matrix<T>& m) const ;
	

private:
	
	/**
	 * @brief   Transform
	 *
	 * @param   m   To transform
	 * @param   bw  Backward: true, Forward: false
	 * @return      Transform
	 */
	Matrix<cxfl> 
	Transform    (const Matrix<cxfl>& m, const bool& bw) const ;

	wlfamily m_wf;                 /**< @brief wavelet family */

	size_t  m_sz;                  /**< @brief data size */
	size_t  m_sl;                  /**< @brief side length */

	double* m_re;                  /**< @brief Real store */
	double* m_im;                  /**< @brief Imag store */
	
	gsl_wavelet_workspace* m_work; /**< @brief Work space */
	gsl_wavelet           *m_w;    /**< @brief Wavelet    */
	
};

#endif
