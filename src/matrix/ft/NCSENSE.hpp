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

#include "NFFT.hpp"
#include "CX.hpp"
#include "SEM.hpp"
#include "MRI.hpp"
#include "Lapack.hpp"
#include "tinyxml.h"

#include "Workspace.hpp"

/**
 * @brief Non-Cartesian SENSE<br/>
 *        According Pruessmann et al. (2001). MRM, 46(4), 638-51.
 *
 */
template <class T>
class NCSENSE : public FT<T> {
	
public:

	/**
	 * @brief         Default constructor
	 */
	NCSENSE() : m_initialised (false),
                m_dim(2),
                m_cgiter(100),
                m_fts(0),
                m_cgeps (1.0e-6),
                m_nc (8),
                m_nk (1024),
                m_nr (4096),
                m_lambda (1.0e-6),
                m_verbose (false) {}
    
    
	/**
	 * @brief          Construct with parameters
	 *
	 * @param  params  Configuration parameters
	 */
	NCSENSE        (const Params& params)
              : FT<T>::FT(params),
                m_nc(0),
                m_nk(0),
                m_nr(0),
                m_fts(0),
                m_cgiter(0),
                m_dim(0),
                m_initialised(false),
                m_cgeps(0.0),
                m_lambda(0.0),
                m_verbose (false) {
        
		T fteps = 7.0e-4, alpha = 1.0;
		size_t ftiter = 3, m = 1;
        
		if (params.exists("fteps"))
			fteps = boost::any_cast<T>(params["fteps"]);
		if (params.exists("alpha"))
			alpha = boost::any_cast<T>(params["alpha"]);
		if (params.exists("ftiter"))
			ftiter = boost::any_cast<size_t>(params["ftiter"]);
		if (params.exists("m"))
			m      = boost::any_cast<size_t>(params["m"]);

		Workspace& ws = Workspace::Instance();
		Matrix<T> b0;

		m_smname = params.Get<std::string>("sens_maps");
		m_wname  = params.Get<std::string>("weights_name");
        m_verbose = params.Get<bool>("verbose");

		if (params.exists("phase_cor"))
			ws.GetMatrix(params.Get<std::string>("phase_cor"), m_pc);
		if (params.exists("b0"))
			ws.GetMatrix(params.Get<std::string>("b0"), b0);

		m_dim = ndims(m_sm)-1;
		Matrix<size_t> ms (m_dim,1);
		for (size_t i = 0; i < m_dim; i++)
			ms[i] = size(m_sm,i);

		m_cgiter = params.Get<size_t>("cgiter");
		m_cgeps  = params.Get<size_t>("cgeps");
		m_lambda = params.Get<size_t>("lambda");
		m_nk     = params.Get<size_t>("nk");

		printf ("  Initialising NCSENSE:\n");
		printf ("  Signal nodes: %li\n", m_nk);
		printf ("  CG: eps(%.3e) iter(%li) lambda(%.3e)\n", m_cgeps, m_cgiter, m_lambda);
		printf ("  FT: eps(%.3e) iter(%li) m(%li) alpha(%.3e)\n", fteps, ftiter, m, alpha);

		int np = 1;

#pragma omp parallel default (shared)
		{
			if (params.exists("np"))
				np = boost::any_cast<int>(params["np"]);
			else
				np = omp_get_num_threads ();
		}

		m_fts = new NFFT<T>* [np];
		for (size_t i = 0; i < np; i++)
			m_fts[i] = new NFFT<T> (ms, m_nk, m, alpha, b0, m_pc, fteps, ftiter);

		m_ic     = IntensityMap (m_sm);
		m_initialised = true;
		printf ("  ...done.\n\n");

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
	NCSENSE (const Matrix< std::complex<T> >& sens, const size_t& nk, const T& cgeps, const size_t& cgiter,
			 const T& lambda = 0.0, const T& fteps = 7.0e-4, const size_t& ftiter = 3,
			 const size_t& m = 1, const T& alpha = 1.0, const Matrix<T>& b0 = Matrix<T>(1),
			 const Matrix< std::complex<T> >& pc = Matrix< std::complex<T> >(1)) :
		m_nc(8),
		m_nk(4096),
		m_nr(1024) {
		

		m_dim = ndims(sens)-1;
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
	virtual ~NCSENSE () {
		
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
	KSpace (const Matrix<T>& k) {
		
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
	Weights (const Matrix<T>& w) {
		
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
	virtual Matrix< std::complex<T> >
	Trafo       (const Matrix< std::complex<T> >& m) const {

		Matrix< std::complex<T> > tmp = m / m_ic;
		return E (tmp, m_sm, m_fts);

	}

	
	/**
	 * @brief    Forward transform
	 *
	 * @param  m     To transform
	 * @param  sens  Sensitivities
	 * @param  recal Recompute intensity connection
	 *
	 * @return   Transform
	 */
	virtual Matrix< std::complex<T> >
	Trafo       (const Matrix< std::complex<T> >& m, const Matrix< std::complex<T> >& sens, const bool& recal = true) const {

		return E (m / ((recal) ? IntensityMap (sens) : m_ic), sens, m_fts);

	}

	
	/**
	 * @brief Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix< std::complex<T> >
	Adjoint     (const Matrix< std::complex<T> >& m) const {

		return this->Adjoint (m, m_sm, false);

	}
	
	
	/**
	 * @brief Backward transform
	 *
	 * @param  m     To transform
	 * @param  sens  Sensitivities
	 * @param  recal Recompute intensity correction (default: true)
	 *
	 * @return   Transform
	 */
	virtual Matrix< std::complex<T> >
	Adjoint (const Matrix< std::complex<T> >& m,
			 const Matrix< std::complex<T> >& sens,
			 const bool recal = true) const {

		std::complex<T> ts;
		T rn, rno, xn;
		Matrix< std::complex<T> > p, r, x, q;
		vector<T> res;
        std::vector< Matrix<cxfl> > vc;

		p = EH (m, sens, m_fts) * ((recal) ? IntensityMap (sens) : m_ic);
		r = p;
		x = zeros< std::complex<T> >(size(p));
		
		xn = pow(norm(p), 2.0);
		rn = xn;
		
		for (size_t i = 0; i < m_cgiter; i++) {
			
			res.push_back(rn/xn);
			
			if (std::isnan(res.at(i)) || res.at(i) <= m_cgeps) 
				break;

 			printf ("    %03lu %.7f\n", i, res.at(i)); fflush (stdout);

			q   = EH (E  (p * ((recal) ? IntensityMap (sens) : m_ic), sens, m_fts), sens, m_fts) *
					((recal) ? IntensityMap (sens) : m_ic);
			if (m_lambda)
				q  += m_lambda * p;
			
			ts  = rn / p.dotc(q);
			x  += ts * p;
            if (m_verbose)
                vc.push_back (x * ((recal) ? IntensityMap (sens) : m_ic));
			r  -= ts * q;
			rno = rn;
			rn  = pow(norm(r), 2.0);
			p  *= rn / rno;
			p  += r;

		}

		printf ("\n");

        if (m_verbose) {
            size_t cpsz = numel(x);
            x = zeros<std::complex<T> > (size(x,0), size(x,1), (m_dim == 3) ? size(x,2) : 1, vc.size());
            for (size_t i = 0; i < vc.size(); i++)
                memcpy (&x[i*cpsz], &(vc[i][0]), cpsz*sizeof(std::complex<T>));
            return x;
        } else        
            return x * ((recal) ? IntensityMap (sens) : m_ic);

	}
	
	
	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix< std::complex<T> >
	operator* (const Matrix< std::complex<T> >& m) const {
		return Trafo(m);
	}
	

	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix< std::complex<T> >
	operator->* (const Matrix< std::complex<T> >& m) const {
		return Adjoint (m);
	}

	
	
private:

	NFFT<T>** m_fts;         /**< Non-Cartesian FT operators (Multi-Core?) */
	bool      m_initialised; /**< All initialised? */
    bool      m_verbose;

	Matrix< std::complex<T> > m_sm;          /**< Sensitivities */
	Matrix<T> m_ic;     /**< Intensity correction I(r) */
	Matrix< std::complex<T> > m_pc; /**< @brief Correction phase */

	std::string m_smname;
	std::string m_wname;

	size_t    m_dim;            /**< Image dimensions {2,3} */
	size_t    m_nr;             /**< # spatial image positions */
	size_t    m_nk;             /**< # K-space points */
	size_t    m_nc;             /**< # Receive channels */

	size_t    m_cgiter;         /**< Max # CG iterations */
	double    m_cgeps;          /**< Convergence limit */
	double    m_lambda;         /**< Tikhonov weight */

};


#endif

