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

#ifndef __NCSENSE_HPP__
#define __NCSENSE_HPP__

#include "nfftstub.h"
#include "NFFT.hpp"
#include "CX.hpp"

/**
 * @brief Matrix templated ND non-equidistand Fourier transform with NCSENSE 3 (TU Chemnitz)
 */
template <class T>
class NCSENSE : public FT<T> {
	
public:

	NCSENSE() : m_initialised (false) {};

	/**
	 * @brief          Construct NCSENSE plans for forward and backward FT with credentials
	 * 
	 * @param  sens    Sensitivity maps if imsize
	 * @param  imsize  Matrix of side length of the image space
	 * @param  nk      # k-space points
	 * @param  m       Spatial cut-off of FT
	 * @param  alpha   Oversampling factor
	 * @param  b0      Off-resonance maps if available
	 * @param  pc      Phase correction applied before forward or after adjoint transforms (default: empty)
	 * @param  eps     Convergence criterium for inverse transform (default: 1.0e-7)
	 * @param  maxit   Maximum # NCSENSE iterations (default: 3)
	 */
	NCSENSE (const Matrix<T> sens, const size_t& nk, const size_t m = 1, 
			 const double alpha = 1.0, const Matrix<double> b0 = Matrix<double>(1), 
			 const Matrix<T> pc = Matrix<T>(1), 
			 const double eps = 7.0e-4, const size_t maxit = 3);
	
	
	/**
	 * @brief        Clean up and destruct
	 */ 
	~NCSENSE () {

	int np = 1;

#pragma omp parallel default (shared)
	{
		np = omp_get_num_threads ();
	}	

	if (m_initialised)
		for (size_t i = 0; i < np; i++)
			delete m_nffts[i];

	}
	
	
	/**
	 * @brief      Assign k-space 
	 * 
	 * @param  k   Kspace trajectory
	 */
	void
	KSpace (const Matrix<double>& k) {

#pragma omp parallel 
		{
			
			for (size_t i = 0; i < omp_get_num_threads (); i++)
				m_nffts[i]->KSpace(k);
			
		}
		
	}
	

	/**
	 * @brief      Assign k-space weigths (jacobian of trajectory with regards to time) 
	 * 
	 * @param  w   Weights
	 */
	void
	Weights (const Matrix<double>& w) {
		
#pragma omp parallel 
		{
			
			for (size_t i = 0; i < omp_get_num_threads (); i++)
				m_nffts[i]->Weights(w);
			
		}
		
	}


	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T> 
	Trafo       (const Matrix<T>& m) const {};
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T> 
	Adjoint     (const Matrix<T>& m) const {};
	
	
private:

	NFFT<T>** m_nffts;
	bool      m_initialised;


};


template<>
NCSENSE<cxdb>::NCSENSE (const Matrix<cxdb> sens, const size_t& nk, const size_t m, 
						const double alpha, const Matrix<double> b0, const Matrix<cxdb> pc, 
						const double eps, const size_t maxit) : m_initialised (false) {

	size_t dim = (size(sens,2) == 1) ? 3 : 2;
	Matrix<size_t> ms (dim,1);
	for (size_t i = 0; i < dim; i++)
		ms[i] = size(sens,i);
	
	int np = 1;

#pragma omp parallel default (shared)
	{
		np = omp_get_num_threads ();
	}	

	m_nffts = new NFFT<cxdb>* [np];
	
	for (size_t i = 0; i < np; i++)
		m_nffts[i] = new NFFT<cxdb> (ms, nk, m, alpha);
	
	m_initialised = true;

}




#endif
