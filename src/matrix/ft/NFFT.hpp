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
    typedef typename std::vector<std::complex<double> >::iterator it;
	
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
		m_rank (0) {};

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
		
		for (size_t i = 0; i < 4; i++) {
			m_N[i] = 1; m_n[i] = 1;
		}
		
		m_rank = numel (imsize);
		
		for (size_t i = 0; i < m_rank; i++) {
			m_N[i]   = imsize[i];
			m_n[i]   = ceil (m_N[i]*alpha);
			m_imgsz *= m_N[i];
		}
		
		m_epsilon = eps;
		m_maxit   = maxit;
		
		NFFTTraits<double>::Init (m_rank, m_N, m_M, m_n, m, m_fplan, m_iplan);
		
        m_y = it((std::complex<double>*)m_iplan.y);
        m_f = it((std::complex<double>*)m_iplan.f_hat_iter);

		if (pc.Size() > 1)
			m_have_pc = true;
		
		m_pc   = pc;
		m_cpc  = conj(pc);
		
		m_initialised = true;
		
	}
	
	
	/**
	 * @brief        Clean up and destruct
	 */ 
	~NFFT () {
		
		if (m_initialised) {
			
			NFFTTraits<double>::Finalize (m_fplan, m_iplan);

		}
		
	}


	/**
	 * @brief      Assign k-space 
	 * 
	 * @param  k   Kspace trajectory
	 */
	inline void 
	KSpace (const Matrix<T>& k) {
		
		if (sizeof(T) == sizeof(double))
			memcpy (m_fplan.x, k.Memory(), m_fplan.M_total * m_rank * sizeof(double));
		else 
			for (size_t i = 0; i < m_fplan.M_total * m_rank; i++)
				m_fplan.x[i] = k[i];
		
	}
	
	
	/**
	 * @brief      Assign k-space weigths (jacobian of trajectory with regards to time) 
	 * 
	 * @param  w   Weights
	 */
	inline void 
	Weights (const Matrix<T>& w) {
		
		if (sizeof(T) == sizeof(double))
			memcpy (m_iplan.w, w.Memory(), m_fplan.M_total * sizeof(double));
		else 
			for (size_t i = 0; i < m_fplan.M_total; i++)
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
			memcpy (m_fplan.f_hat, m.Memory(), m_imgsz * sizeof (std::complex<T>));
		else 
			for (size_t i = 0; i < m_imgsz; i++) {
				m_fplan.f_hat[i][0] = creal(m[i]);
				m_fplan.f_hat[i][1] = cimag(m[i]);
			}
		
		NFFTTraits<double>::Trafo (m_fplan);
		
		if (sizeof(T) == sizeof(double))
			memcpy (&out[0], m_fplan.f, m_M * sizeof (std::complex<T>));
		else 
			for (size_t i = 0; i < m_M; i++)
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

        Matrix< std::complex<T> > out (m_N[0], m_N[1], m_N[2], m_N[3]);
		
		if (sizeof(T) == sizeof(double))
			memcpy (m_iplan.y, m.Memory(), m_M * sizeof ( std::complex<T> ));
		else 
			for (size_t i = 0; i < m_M; i++) {
				m_iplan.y[i][0] = creal(m[i]);
				m_iplan.y[i][1] = cimag(m[i]);
			}

		NFFTTraits<double>::ITrafo ((Plan&) m_fplan, (Solver&) m_iplan, m_maxit, m_epsilon);

		if (sizeof(T) == sizeof(double))
			memcpy (&out[0], m_iplan.f_hat_iter, m_imgsz*sizeof ( std::complex<T> ));
		else 
			for (size_t i = 0; i < m_imgsz; i++)
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
	
	bool      m_initialised;      /**< @brief Memory allocated / Plans, well, planned! :)*/
	bool      m_have_pc;

	size_t    m_rank;

	Matrix< std::complex<T> > m_pc;  /**< @brief Phase correction (applied after inverse trafo)*/
	Matrix< std::complex<T> > m_cpc; /**< @brief Phase correction (applied before forward trafo)*/
	
	int       m_N[4];             /**< @brief Image matrix side length (incl. k_{\\omega})*/
	int       m_n[4];             /**< @brief Oversampling */
	int       m_M;                /**< @brief Number of k-space knots */
	int       m_maxit;            /**< @brief Number of Recon iterations (NFFT 3) */
	T         m_epsilon;          /**< @brief Convergence criterium */
	size_t    m_imgsz;
	
	Plan      m_fplan;            /**< nfft  plan */
	Solver    m_iplan;            /**< infft plan */
	
	bool      m_prepared;

    it        m_y;
    it        m_f;


};


#endif




