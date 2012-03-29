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

#include "NFFTStub.hpp"
#include "FT.hpp"
#include "CX.hpp"

/**
 * @brief Matrix templated ND non-equidistand Fourier transform with NFFT 3 (TU Chemnitz)
 */
template <class T>
class NFFT : public FT<T> {
	
public:

	NFFT() : m_initialised (false) {};

	/**
	 * @brief        Construct NFFT plans for forward and backward FT with credentials
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
				 const double alpha = 1.0, const Matrix<T> b0 = Matrix<T>(1), 
				 const Matrix<T> pc = Matrix<T>(1), const double eps = 7.0e-4, 
				 const size_t maxit = 3);
	
	
	/**
	 * @brief        Clean up and destruct
	 */ 
	~NFFT () {
		
		if (m_initialised) {
			
			solver_finalize_complex (&m_iplan);
			nfft_finalize (&m_fplan);

		}
		
	}


	/**
	 * @brief      Assign k-space 
	 * 
	 * @param  k   Kspace trajectory
	 */
	void 
	KSpace (const Matrix<double>& k) {
		
		memcpy (m_fplan.x, k.Data(), numel(k) * sizeof(double));
		
	}
	
	
	/**
	 * @brief      Assign k-space weigths (jacobian of trajectory with regards to time) 
	 * 
	 * @param  w   Weights
	 */
	void 
	Weights (const Matrix<double>& w)  {
		
		memcpy (m_iplan.w, w.Data(), numel(w) * sizeof(double));
		
		nnfft::weights (m_fplan, m_iplan);
		nnfft::psi     (m_fplan);
		
	}
	

	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T> 
	Trafo       (const Matrix<T>& m) const ;
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T> 
	Adjoint     (const Matrix<T>& m) const ;
	
	
private:
	
	bool      m_initialised;      /**< @brief Memory allocated / Plans, well, planned! :)*/
	bool      m_have_pc;

	Matrix<T> m_pc;               /**< @brief Phase correction (applied after inverse and before forward trafos) (double precision)*/
	Matrix<T> m_cpc;              /**< @brief Phase correction (applied after inverse and before forward trafos) (double precision)*/
	
	int       m_N[4];             /**< @brief Image matrix side length (incl. k_\omega)*/
	int       m_n[4];             /**< @brief Oversampling */
	int       m_M;                /**< @brief Number of k-space knots */
	int       m_maxit;            /**< @brief Number of Recon iterations (NFFT 3) */
	double    m_epsilon;          /**< @brief Convergence criterium */
	size_t    m_imgsz;
	
	nfft_plan m_fplan;            /**< nfft  plan */
	solver_plan_complex m_iplan;  /**< infft plan */
	
	bool      m_prepared;

};


template<>
NFFT<cxdb>::NFFT (const Matrix<size_t>& imsize, const size_t& nk, const size_t m, const double alpha, 
				  const Matrix<cxdb> b0, const Matrix<cxdb> pc, const double eps, const size_t maxit)  
  : m_initialised (false) {

	m_M = nk;
	m_imgsz = 1;

	for (size_t i = 0; i < 4; i++) {
		m_N[i] = 1; m_n[i] = 1;
	}
	
	int  rank = numel (imsize);
	
	for (size_t i = 0; i < rank; i++) {
		m_N[i] = imsize[i];
		m_n[i] = ceil (m_N[i]*alpha);
		m_imgsz *= m_N[i];
	}
	
	m_epsilon = eps;
	m_maxit   = maxit;
	
	nnfft::init (rank, m_N, m_M, m_n, m, m_fplan, m_iplan);

	if (pc.Size() > 1)
		m_have_pc = true;

	m_pc   = pc;
	m_cpc  = conj(pc);

	m_initialised = true;

}


template<> inline Matrix<cxdb>
NFFT<cxdb>::Adjoint (const Matrix<cxdb>& in) const {

	Matrix<cxdb> out (m_N[0], m_N[1], m_N[2], m_N[3]);
	
	memcpy (m_iplan.y, in.Data(), m_M * sizeof (cxdb));
	
	nnfft::ift ((nfft_plan&) m_fplan, (solver_plan_complex&) m_iplan, m_maxit, m_epsilon);
	
	memcpy (&out[0], m_iplan.f_hat_iter, m_imgsz*sizeof (cxdb));
	
	return out;

}


template<> inline Matrix<cxdb>
NFFT<cxdb>::Trafo (const Matrix<cxdb>& in) const {

	Matrix<cxdb> out (m_M,1);

	memcpy (m_fplan.f_hat, in.Data(), m_imgsz * sizeof (cxdb));

	nnfft::ft (m_fplan);

	memcpy (&out[0], m_fplan.f, m_M * sizeof (cxdb));
			
	return out;

}


template<> inline Matrix<cxfl>
NFFT<cxfl>::Adjoint (const Matrix<cxfl>& in) const {

	Matrix<cxfl> out (m_N[0],m_N[1],m_N[2],m_N[3]);
	size_t m = m_M-1, i = m_imgsz-1;  

	while (m--) {
		m_iplan.y[m][0] = in[m].real();
		m_iplan.y[m][1] = in[m].imag();
	}

	nnfft::ift ((nfft_plan&) m_fplan, (solver_plan_complex&) m_iplan, m_maxit, m_epsilon);

	while (i--)
		out [i] = cxfl (m_iplan.f_hat_iter[i][0], m_iplan.f_hat_iter[i][1]);
	
	return out;

}


template<> inline Matrix<cxfl>
NFFT<cxfl>::Trafo (const Matrix<cxfl>& in) const {

	Matrix<cxfl> out (m_M,1);

	for (size_t i = 0; i < m_M; i++) {
		m_fplan.f_hat[i][0] = in[i].real();
		m_fplan.f_hat[i][1] = in[i].imag();
	}

	nnfft::ft (m_fplan);

	for (size_t i = 0; i < m_imgsz; i++) 
		out [i] = cxfl (m_fplan.f[i][0], m_fplan.f[i][1]);
			
	return out;

}


	
#endif




