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
#include "MRI.hpp"

/**
 * @brief   Convenience structure for construction
 */
template<class T>
struct NCSParams {
	
	Matrix<T> sens;
	Matrix<T> pc;

	Matrix<double> b0; 

	double cgeps;
	double lambda;
	double fteps; 
	double alpha; 

	size_t nk;
	size_t cgiter; 
	size_t ftiter; 
	size_t m; 

	NCSParams () {
		
		cgeps  = 1.0e-6;
		lambda = 2.0e-6;
		fteps  = 7.0e-4;
		alpha  = 1.0;

		ftiter = 3;
		m      = 1;
		cgiter = 20;
		nk     = 0;
		
	}

};
		

/**
 * @brief Non-Cartesian SENSE<br/>
 *        According Pruessmann et al. (2001). MRM, 46(4), 638-51.
 */
template <class T>
class NCSENSE : public FT<T> {
	
public:

	NCSENSE() : m_initialised (false) {};

	NCSENSE (const NCSParams<T>& ncsp) {

	}


	/**
	 * @brief          Construct NCSENSE plans for forward and backward transform with credentials
	 * 
	 * @param  sens    Sensitivity maps if imsize
	 * @param  nk      # k-space points
	 * @param  cgeps   Convergence limit of descent
	 * @param  cgiter  Maximum # CG iterations
	 * @param  lambda  Tikhonov regularisation (default 0.0)
	 * @param  fteps   NFFT convergence criterium (default 7.0e-4)
	 * @param  ftiter  Maximum # of NFFT gridding iterations (default 3)
	 * @param  m       Spatial cut-off of FT (default 1)
	 * @param  alpha   Oversampling factor (default 1.0)
	 * @param  b0      Off-resonance maps if available (default empty)
	 * @param  pc      Phase correction applied before forward or after adjoint transforms (default: empty)
	 */
	NCSENSE (const Matrix<T>& sens, const size_t& nk, const double& cgeps, const size_t& cgiter, 
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
		m_ic     = IntensityMap (m_sm);
		
		m_initialised = true;
		
		printf ("  ...done.\n\n");

	}
	
	/**
	 * @brief        Clean up and destruct NFFT plans
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
	 * @brief      Assign k-space trajectory
	 * 
	 * @param  k   K-space trajectory
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
	 * @brief      Assign k-space weigths (jacobian of k in t) 
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
		return E (tmp, m_sm, m_fts);

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

		p = EH (m, m_sm, m_fts) * m_ic;
		r = p;
		x = zeros<T>(size(p));
		
		rn = 0.0;
		xn = pow(creal(norm(p)), 2);
		
		for (size_t i = 0; i < m_cgiter; i++) {
			
			rn  = pow(creal(norm(r)), 2);
			res.push_back(rn/xn);
			
			if (std::isnan(res.at(i)) || res.at(i) <= m_cgeps) 
				break;
			
			printf ("    %03lu %.7f\n", i, res.at(i)); fflush (stdout);
			
			q   = E  (p * m_ic, m_sm, m_fts);
			q   = EH (q, m_sm, m_fts) * m_ic;
			
			if (m_lambda)
				q  += m_lambda * p;
			
			ts  = rn;
			ts /= p.dotc(q);
			x  += ts * p;
			r  -= ts * q;
			p  *= pow(creal(norm(r)), 2)/rn;
			p  += r;
			
		}

		printf ("\n");
		return x * m_ic;

	}
	
	
	
	
private:

	NFFT<T>** m_fts;         /**< Non-Cartesian FT operators (Multi-Core?) */
	bool      m_initialised; /**< All initialised? */

	Matrix<T> m_sm;          /**< Sensitivities */
	Matrix<double> m_ic;     /**< Intensity correction I(r) */

	size_t m_dim;            /**< Image dimensions {2,3} */
	size_t m_nr;             /**< # spatial image positions */
	size_t m_nk;             /**< # K-space points */
	size_t m_nc;             /**< # Receive channels */

	size_t m_cgiter;         /**< Max # CG iterations */
	double m_cgeps;          /**< Convergence limit */
	double m_lambda;         /**< Tikhonov weight */

};


#endif

