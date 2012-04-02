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


template <class T> inline static Matrix<double>
IntensityCorrection (const Matrix<T>& sens) {

	size_t dim = ndims(sens);
	size_t nc  = size(sens,dim);
	size_t nr  = numel(sens)/nc;
	
	Matrix<size_t> dims = size (sens);
	dims [dim] = 1;

	Matrix<double> res = zeros<double> (dims);
	
#pragma omp parallel default (shared)
	{		
		
#pragma omp for schedule (guided)
		for (size_t i = 0; i < nr; i++) {
			
			for (size_t j = 0; j < nc; j++)
				res[i] += (sens(i+j*nr) * conj(sens(i+j*nr))).real();
			
			res[i] = 1.0 / (sqrt (res[i]) + 1.0e-10);
			
		}
		
	}

	return res;

}

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
			 const double& lambda = 0.0, const double& fteps = 7.0e-4, const size_t& ftiter = 3, 
			 const size_t& m = 1, const double& alpha = 1.0, const Matrix<double>& b0 = Matrix<double>(1), 
			 const Matrix<T>& pc = Matrix<T>(1)) {
		

		m_dim = ndims(sens);
		Matrix<size_t> ms (m_dim,1);
		for (size_t i = 0; i < m_dim; i++)
			ms[i] = size(sens,i);
		
		printf ("  Initialising NCSENSE:\n");
		printf ("  Signal nodes: %li\n", nk);
		printf ("  CG: eps(%.3e) iter(%li) lambda(%.3e)\n", cgeps, cgiter, lambda);
		printf ("  FT: eps(%.3e) iter(%li) m(%li) alpha(%.3e)\n", fteps, ftiter, m, alpha);

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

		m_sm     = sens;
		m_ic     = IntensityCorrection (m_sm);
		
		m_initialised = true;
		
		printf ("  ...done.\n\n");

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

		Matrix<T> tmp = m * m_ic;
		return E (tmp, m_sm, m_fts, m_dim);

	}

	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	Matrix<T> 
	Adjoint     (const Matrix<T>& m) const {

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
			
			if (std::isnan(res.at(i)) || res.at(i) <= m_cgeps) 
				break;
			
			printf ("    %03lu %.7f\n", i, res.at(i)); fflush (stdout);
			
			q   = E  (p * m_ic, m_sm, m_fts, m_dim);
			q   = EH (q, m_sm, m_fts, m_dim) * m_ic;
			
			if (m_lambda)
				q  += m_lambda * p;
			
			ts  = rn;
			ts /= p.dotc(q);
			x  += ts * p;
			r  -= ts * q;
			p  *= pow(creal(Norm(r)), 2)/rn;
			p  += r;
			
		}

		printf ("\n");
		return x * m_ic;

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

/*
	printf ("  intialising nfft::init (%i, {%i, %i, %i}, %i, {%i, %i, %i}, %i, *, *, %.9f)\n", 
			m_dim, 
			m_N[0], m_N[1], m_N[2],
			m_M,
			m_n[0], m_n[1], m_n[2],
			m,
			m_fteps);

	for (int i = 0; i < NTHREADS || i < m_Nc; i++)
		nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan[i], &m_iplan[i]);
	// --------------------------------------

	*/
