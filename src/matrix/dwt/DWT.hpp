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


# include "Matrix.hpp"
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
    DWT (const size_t& sl, const wlfamily& wf = WL_DAUBECHIES, const size_t& wm = 4) {
        
        // Checks missing !!!
        
        m_wf = wf;
        
        if (m_wf > ID) {
            
            m_sl   = sl;
            m_sz   = sl*sl;
            
            switch (wf) {
                case WL_DAUBECHIES         : m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies,          wm); break;
                case WL_DAUBECHIES_CENTERED: m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies_centered, wm); break;
                case WL_HAAR               : m_w = gsl_wavelet_alloc (gsl_wavelet_haar,                wm); break;
                case WL_HAAR_CENTERED      : m_w = gsl_wavelet_alloc (gsl_wavelet_haar_centered,       wm); break;
                case WL_BSPLINE            : m_w = gsl_wavelet_alloc (gsl_wavelet_bspline,             wm); break;
                case WL_BSPLINE_CENTERED   : m_w = gsl_wavelet_alloc (gsl_wavelet_bspline_centered,    wm); break;
                default                    : m_w = gsl_wavelet_alloc (gsl_wavelet_daubechies,          wm); break;
            }
            
            m_work = gsl_wavelet_workspace_alloc (m_sz);
            m_iwork = gsl_wavelet_workspace_alloc (m_sz);

            m_re = VECTOR_CONSTR(Type,m_sz);
            if (typeid(T) == cxfl_type || typeid(T) == cxdb_type)
                m_im = VECTOR_CONSTR(Type,m_sz);
            
        }
        
    }
    

    virtual
    ~DWT () {
        
        if (m_wf != ID) {
            gsl_wavelet_workspace_free (m_work);
            gsl_wavelet_workspace_free (m_iwork);
        }

        
    }
    
    
	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
    inline Matrix<T>
	Trafo        (const Matrix<T>& m) {
		return Transform (m);
	}
	

	/**
	 * @brief    Adjoint transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline Matrix<T>
	Adjoint      (const Matrix<T>& m) {
		return Transform (m, true);
	}
	

	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
    inline Matrix<T>
	operator*    (const Matrix<T>& m) {
		return Trafo(m);
	}
	

	/**
	 * @brief    Adjoint transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
    inline Matrix<T>
	operator->* (const Matrix<T>& m) {
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
	inline Matrix<T>
	Transform    (const Matrix<T>& m, const bool& bw = false)  {
        
        Matrix<T> res = m;
        
        if (m_wf > ID) {
            
            // Checks missing !!!
            
            DWTTraits<T>::prepare (res, m_re, m_im, m_sl);
            
            omp_set_num_threads (2);

            if (bw) {
                
#pragma omp parallel default (shared)
            	{
#pragma omp sections nowait
            		{
#pragma omp section
            			if (!(DWTTraits<T>::nstransform_inverse (m_w, m_re, m_sl, m_sl, m_sl, m_work) == GSL_SUCCESS))
            				printf ("Wavelet transform for real part failed\n.");
#pragma omp section
            			if (!(DWTTraits<T>::nstransform_inverse (m_w, m_im, m_sl, m_sl, m_sl, m_iwork, true) == GSL_SUCCESS))
            				printf ("Wavelet transform for imaginary part failed\n.");
            		}
            	}

            } else {
                
#pragma omp parallel default (shared)
            	{
#pragma omp sections nowait
            		{
#pragma omp section
						if (!(DWTTraits<T>::nstransform_forward (m_w, m_re, m_sl, m_sl, m_sl, m_work) == GSL_SUCCESS))
							printf ("Wavelet transform for real part failed\n.");
#pragma omp section
						if (!(DWTTraits<T>::nstransform_forward (m_w, m_im, m_sl, m_sl, m_sl, m_iwork, true) == GSL_SUCCESS))
                			printf ("Wavelet transform for imaginary part failed\n.");
                	}
                }

            }
            
            DWTTraits<T>::finalize (res, m_re, m_im, m_sl);
            
        }
        
        return res;
        
    }
    
	wlfamily m_wf;                 /**< @brief wavelet family */
    
	size_t  m_sz;                  /**< @brief data size */
	size_t  m_sl;                  /**< @brief side length */
    
	VECTOR_TYPE(double) m_re;      /**< @brief Real store */
	VECTOR_TYPE(double) m_im;      /**< @brief Imag store */
	
	gsl_wavelet_workspace *m_work; /**< @brief Work space */
	gsl_wavelet_workspace *m_iwork;
	gsl_wavelet           *m_w;    /**< @brief Wavelet    */
    
	
};


# endif // __DWT_HPP__
