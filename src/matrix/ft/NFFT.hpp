/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                                  Forschungszentrum Juelich, Germany
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

#ifndef __NFFT_HPP__
#define __NFFT_HPP__

#include "NFFTTraits.hpp"
#include "Algos.hpp"
#include "FT.hpp"
#include "CX.hpp"
#include "Creators.hpp"
/**
 * @brief Matrix templated ND non-equidistand Fourier transform with NFFT 3 (TU Chemnitz)<br/>
 *        Double and single precision
 */
template <class T>
class NFFT : public FT<T> {


	typedef typename NFFTTraits<double>::Plan   Plan;
	typedef typename NFFTTraits<double>::Solver Solver;
  //typedef typename std::vector<std::complex<double> >::iterator it;
    typedef typename std::complex<T> CT;
	
public:

	/**
	 * @brief         Default constructor
	 */
	NFFT() :
		m_initialised (false),
		m_have_pc (false),
		m_imgsz (0),
		m_M (0),
		m_maxit (0),
		m_prepared (false),
		m_rank (0),
		m_m(0) {};

	/**
	 * @brief          Construct NFFT plans for forward and backward FT with credentials
	 * 
	 * @param  imsize  Matrix of side length of the image space
	 * @param  nk      # k-space points
	 * @param  m       Spatial cut-off of FT
	 * @param  alpha   Oversampling factor
	 * @param  b0      Off-resonance maps if available
	 * @param  pc      Phase correction applied before forward or after adjoint transforms (default: empty)
	 * @param  eps     Convergence criterium for inverse transform (default: 1.0e-7)
	 * @param  maxit   Maximum # NFFT iterations (default: 3)
	 */
	NFFT        (const Matrix<size_t>& imsize, const size_t& nk, const size_t m = 1, 
				 const T alpha = 1.0, const Matrix<T> b0 = Matrix<T>(1),
				 const Matrix< std::complex<T> > pc = Matrix< std::complex<T> >(1),
				 const T eps = 7.0e-4, const size_t maxit = 1) : m_prepared (false) {
		
		m_M     = nk;
		m_imgsz = 1;
		m_m = m;
		
		m_N = imsize.Container();
		m_n = m_N;//ceil (alpha*m_N);

		m_rank = numel(imsize);
		
		m_imgsz = prod(m_N);

		m_epsilon = eps;
		m_maxit   = maxit;
		
		NFFTTraits<double>::Init (m_N, m_M, m_n, m_m, m_fplan, m_iplan);
		
        m_y = m_iplan.y;
        m_f = m_iplan.f_hat_iter;

		if (pc.Size() > 1)
			m_have_pc = true;
		
		m_pc   = pc;
		m_cpc  = conj(pc);
		
		m_initialised = true;
		
	}
	
	
	/**
	 * @brief Copy conctructor
	 */
	NFFT (const NFFT<T>& ft) {
		*this = ft;
	}


	/**
	 * @brief        Clean up and destruct
	 */ 
	~NFFT () {
		
		if (m_initialised)
			NFFTTraits<double>::Finalize (m_fplan, m_iplan);

	}


	/**
	 * @brief 	Assignement
	 */
	inline NFFT<T>& operator= (const NFFT<T>& ft) {

		m_initialised = ft.m_initialised;
		m_have_pc     = ft.m_have_pc;
		m_rank        = ft.m_rank;
		m_pc          = ft.m_pc;
		m_cpc         = ft.m_cpc;
		m_N           = ft.m_N;
		m_n           = ft.m_n;
		m_M           = ft.m_M;
		m_maxit       = ft.m_maxit;
		m_epsilon     = ft.m_epsilon;
		m_imgsz       = ft.m_imgsz;
		m_m           = ft.m_m;
		NFFTTraits<double>::Init (m_N, m_M, m_n, m_m, m_fplan, m_iplan);
		m_prepared    = ft.m_prepared;
	    m_y           = ft.m_y;
	    m_f           = ft.m_f;
	    return *this;
		
	}

	/**
	 * @brief      Assign k-space 
	 * 
	 * @param  k   Kspace trajectory
	 */
	inline void 
	KSpace (const Matrix<T>& k) {
		
		size_t cpsz = k.Size();
		assert (cpsz = m_fplan.M_total * m_rank);

		if (sizeof(T) == sizeof(double))
			memcpy (m_fplan.x, k.Ptr(), cpsz * sizeof(double));
		else 
			for (size_t i = 0; i < cpsz; ++i)
				m_fplan.x[i] = k[i];
		
	}
	
	
	/**
	 * @brief      Assign k-space weigths (jacobian of trajectory with regards to time) 
	 * 
	 * @param  w   Weights
	 */
	inline void 
	Weights (const Matrix<T>& w) {
		
		size_t cpsz = w.Size();
		assert (cpsz = m_fplan.M_total);

		if (sizeof(T) == sizeof(double))
			memcpy (m_iplan.w, w.Ptr(), cpsz * sizeof(double));
		else 
			for (size_t i = 0; i < cpsz; ++i)
				m_iplan.w[i] = w[i];
		
		NFFTTraits<double>::Weights (m_fplan, m_iplan);
		NFFTTraits<double>::Psi     (m_fplan);
		
	}
	
	
	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix< std::complex<T> >
	Trafo       (const Matrix< std::complex<T> >& m) const {

		Matrix< std::complex<T> > out (m_M,1);
		
		if (sizeof(T) == sizeof(double))
			memcpy (m_fplan.f_hat, m.Ptr(), m_imgsz * sizeof (std::complex<T>));
		else 
			for (size_t i = 0; i < m_imgsz; ++i) {
				m_fplan.f_hat[i][0] = real(m[i]);
				m_fplan.f_hat[i][1] = imag(m[i]);
			}
		
		NFFTTraits<double>::Trafo (m_fplan);
		
		if (sizeof(T) == sizeof(double))
			memcpy (&out[0], m_fplan.f, m_M * sizeof (std::complex<T>));
		else 
			for (size_t i = 0; i < m_M; ++i)
				out[i] = std::complex<T>(m_fplan.f[i][0],m_fplan.f[i][1]);
		
		return out;

	}
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix< std::complex<T> >
	Adjoint     (const Matrix< std::complex<T> >& m) const {

        Matrix< std::complex<T> > out (m_N);
		
		if (sizeof(T) == sizeof(double))
			memcpy (m_iplan.y, m.Ptr(), m_M * sizeof ( std::complex<T> ));
		else 
			for (size_t i = 0; i < m_M; ++i) {
				m_iplan.y[i][0] = real(m[i]);
				m_iplan.y[i][1] = imag(m[i]);
			}

		NFFTTraits<double>::ITrafo ((Plan&) m_fplan, (Solver&) m_iplan, m_maxit, m_epsilon);

		if (sizeof(T) == sizeof(double))
			memcpy (&out[0], m_iplan.f_hat_iter, m_imgsz*sizeof ( std::complex<T> ));
		else 
			for (size_t i = 0; i < m_imgsz; ++i)
				out[i] = std::complex<T>(m_iplan.f_hat_iter[i][0],m_iplan.f_hat_iter[i][1]);
		
		return out;
		
	}


	/**
	 * @brief     NFFT plan
	 *
	 * @return    Plan
	 */
	Plan*
	FPlan         () {

		return &m_fplan;

	}
	
	/**
	 * @brief     NFFT plan
	 *
	 * @return    Plan
	 */
	Solver*
	IPlan         () {
		
		return &m_iplan;
		
	}
	
private:
	
	bool       m_initialised;   /**< @brief Memory allocated / Plans, well, planned! :)*/
	bool       m_have_pc;

	size_t     m_rank;

	Matrix<CT> m_pc;            /**< @brief Phase correction (applied after inverse trafo)*/
	Matrix<CT> m_cpc;           /**< @brief Phase correction (applied before forward trafo)*/
	
	container<size_t> m_N;      /**< @brief Image matrix side length (incl. k_{\\omega})*/
	container<size_t> m_n;      /**< @brief Oversampling */

	size_t     m_M;             /**< @brief Number of k-space knots */
	size_t     m_maxit;         /**< @brief Number of Recon iterations (NFFT 3) */
	T          m_epsilon;       /**< @brief Convergence criterium */
	size_t     m_imgsz;
	
	Plan       m_fplan;         /**< nfft  plan */
	Solver     m_iplan;         /**< infft plan */
	
	bool       m_prepared;

    fftw_complex *         m_y;
    fftw_complex *         m_f;

    size_t     m_m;

};


#endif




