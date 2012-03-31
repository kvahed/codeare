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
#include "SEM.hpp"

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
	NCSENSE (const Matrix<T> sens, const size_t& nk, const double& cgeps, const size_t& cgiter, 
			 const double& lambda = 0.0, const size_t& m = 1, const double& alpha = 1.0, 
			 const Matrix<double>& b0 = Matrix<double>(1), const Matrix<T>& pc = Matrix<T>(1), 
			 const double& fteps = 7.0e-4, const size_t& ftiter = 1) {
		
		size_t dim = (size(sens,2) == 1) ? 3 : 2;
		Matrix<size_t> ms (dim,1);
		for (size_t i = 0; i < dim; i++)
			ms[i] = size(sens,i);
		
		int np = 1;
		
#pragma omp parallel default (shared)
		{
			np = omp_get_num_threads ();
		}	
		
		m_fts = new NFFT<T>* [np];
		
		for (size_t i = 0; i < np; i++)
			m_fts[i] = new NFFT<T> (ms, nk, m, alpha, b0, pc, fteps, ftiter);
		
		m_cgiter = cgiter;
		m_cgeps  = cgeps;
		m_lambda = lambda;
		
		printf ("Initialised NC SENSE: %i %e %e\n", m_cgiter, m_cgeps, m_lambda);
		
		m_initialised = true;
		
	}
	
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
			delete m_fts[i];

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
				m_fts[i]->KSpace(k);
			
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
				m_fts[i]->Weights(w);
			
		}
		
	}


	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T> 
	Trafo       (const Matrix<T>& m) const {

		T ts;
		double rn, xn;
		Matrix<T> p, r, x, q;
		vector<double> res;
		
		p = EH (m, m_sm, m_fts, m_dim) * m_ic;
		r = p;
		x = zeros<T>(size(p));
		
		rn = 0.0;
		xn = pow(creal(Norm(p)), 2);
		
		for (size_t i = 0; i < m_cgiter; i++) {
			
			rn  = pow(creal(Norm(r)), 2);
			res.push_back(rn/xn);
			
			if (std::isnan(res.at(i)) || res.at(i) <= m_cgeps) break;
			
			if (i % 5 == 0 && i > 0) printf ("\n");
			printf ("    %03lu %.7f", i, res.at(i));
			
			p  *= m_ic;
			q   = E  (p, m_sm, m_fts, m_dim);
			q   = EH (q, m_sm, m_fts, m_dim);
			q  *= m_ic;
			
			if (m_lambda)
				q  += m_lambda * p;
			
			ts  = rn;
			ts /= p.dotc(q);
			x  += (p * ts);
			r  -= (q * ts);
			p  *= pow(creal(Norm(r)), 2)/rn;
			p  += r;
			
		}
		
		return x * m_ic;

	}
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T> 
	Adjoint     (const Matrix<T>& m) const {

		Matrix<T> tmp = m * m_ic;
		return E (tmp, m_sm, m_fts, m_dim);

	}
	
	
	
private:

	NFFT<T>** m_fts;
	bool      m_initialised;

	Matrix<T> m_sm;
	Matrix<double> m_ic;

	size_t m_dim;
	size_t m_nr;
	size_t m_nk;
	size_t m_nc;

	size_t m_cgiter;
	double m_cgeps;
	double m_lambda;

};


#endif
