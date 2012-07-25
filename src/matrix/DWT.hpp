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


# ifndef __DWT_HPP__

# define __DWT_HPP__


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


/********************
 ** matrix headers **
 ********************/
# include "Matrix.hpp"

/****************
 ** DWT traits **
 ****************/
# include "DWTTraits.hpp"



/**
 * @brief 2D Discrete wavelet transform for Matrix template (from GSL)
 */
template <class T>
class DWT {


	typedef typename DWTTraits<T>::Type Type;


public:


	/**
	 * @brief Construct 2D Wavelet transform with wavelet class and side length
	 *
	 * @param  sl      Side length
	 * @param  wf      Wavelet family (default none, i.e. ID)
	 * @param  wm      Familty member (default 4)
	 */
	DWT (const size_t& sl, const wlfamily& wf = ID, const size_t& wm = 4);
	
	
	virtual 
	~DWT();


	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
  Matrix<T>
	Trafo        (const Matrix<T>& m) const {

		return Transform (m);

	}
	

	/**
	 * @brief    Adjoint transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T>
	Adjoint      (const Matrix<T>& m) const {

		return Transform (m, true);

	}
	

	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
  Matrix<T>
	operator*    (const Matrix<T>& m) const {

		return Trafo(m);

	}
	

	/**
	 * @brief    Adjoint transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
  Matrix<T>
	operator->* (const Matrix<T>& m) const {

		return Adjoint(m);

	}
	
	
private:
	
	/**
	 * @brief   Transform
	 *
	 * @param   m   To transform
	 * @param   bw  Backward: true, Forward: false
	 * @return      Transform
	 */
	Matrix<T>
	Transform    (const Matrix<T>& m, const bool& bw = false) const ;

	wlfamily m_wf;                 /**< @brief wavelet family */

	size_t  m_sz;                  /**< @brief data size */
	size_t  m_sl;                  /**< @brief side length */

	Type * m_re;                   /**< @brief Real store */
	Type * m_im;                   /**< @brief Imag store */
	
	gsl_wavelet_workspace *m_work; /**< @brief Work space */
	gsl_wavelet           *m_w;    /**< @brief Wavelet    */

	
};



/*****************
 ** definitions **
 *****************/


template <class T>
DWT<T>::DWT (const size_t& sl, const wlfamily& wf, const size_t& wm) {

	// Checks missing !!!

	m_wf = wf;

	if (m_wf > ID) {

		m_sl   = sl;
		m_sz   = sl*sl;

		switch (wf) 
			{
			case WL_DAUBECHIES         : m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies,          wm); break;
			case WL_DAUBECHIES_CENTERED: m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies_centered, wm); break;
			case WL_HAAR               : m_w = gsl_wavelet_alloc (gsl_wavelet_haar,                wm); break;
			case WL_HAAR_CENTERED      : m_w = gsl_wavelet_alloc (gsl_wavelet_haar_centered,       wm); break;
			case WL_BSPLINE            : m_w = gsl_wavelet_alloc (gsl_wavelet_bspline,             wm); break;
			case WL_BSPLINE_CENTERED   : m_w = gsl_wavelet_alloc (gsl_wavelet_bspline_centered,    wm); break;
			default                    : m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies,          wm); break;
			}

		m_work = gsl_wavelet_workspace_alloc (m_sz);
		
		m_re = DWTTraits<T>::Malloc (m_sz);
		m_im = DWTTraits<T>::Malloc (m_sz);

	}

}


template <class T>
DWT<T>::~DWT () {

	if (m_wf != ID) {

		free (m_re);
		free (m_im);
		
		gsl_wavelet_workspace_free (m_work);

	}

}


template <class T>
Matrix<T>
DWT<T>::Transform (const Matrix<T>& m, const bool& bw) const {
	
	Matrix<T> res = m;
	
	if (m_wf > ID) {
		
		// Checks missing !!!
		
		DWTTraits<T>::Prepare (res, m_re, m_im, m_sl);
		
		if (bw) {
			
			if (!(DWTTraits<T>::NSTransformInverse (m_w, m_re, m_sl, m_sl, m_sl, m_work) == GSL_SUCCESS))
				printf ("Wavelet transform for real part failed\n.");
			if (!(DWTTraits<T>::NSTransformInverse (m_w, m_im, m_sl, m_sl, m_sl, m_work, true) == GSL_SUCCESS))
				printf ("Wavelet transform for imaginary part failed\n.");
			
		} else {
			
			if (!(DWTTraits<T>::NSTransformForward (m_w, m_re, m_sl, m_sl, m_sl, m_work) == GSL_SUCCESS))
				printf ("Wavelet transform for real part failed\n.");
			if (!(DWTTraits<T>::NSTransformForward (m_w, m_im, m_sl, m_sl, m_sl, m_work, true) == GSL_SUCCESS))
				printf ("Wavelet transform for imaginary part failed\n.");
				
		}
		
		DWTTraits<T>::Finalize (res, m_re, m_im, m_sl);
		
	}
	
	return res;
	
}



# endif // __DWT_HPP__
