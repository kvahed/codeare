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
#include "IOContext.hpp"

#include "Workspace.hpp"

#include <numeric>

/**
 * @brief Non-Cartesian SENSE<br/>
 *        According Pruessmann et al. (2001). MRM, 46(4), 638-51.
 *
 */
template <class T>
class NCSENSE : public FT<T> {

    typedef std::complex<T> CT;
	
public:

	/**
	 * @brief         Default constructor
	 */
	NCSENSE() : m_initialised (false),
                m_cgiter(30),
                m_cgeps (1.0e-6),
                m_lambda (1.0e-6),
                m_verbose (false),
                m_np(0) {}
    
    
	/**
	 * @brief          Construct with parameters
	 *
	 * @param  params  Configuration parameters
	 */
	NCSENSE        (const Params& params)
              : FT<T>::FT(params),
                m_cgiter(30),
                m_initialised(false),
                m_cgeps(1.0e-6),
                m_lambda(1.0e-6),
                m_verbose (false),
		m_np(0) {

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

		assert (params.exists("sens_maps"));
		m_smname = params.Get<std::string>("sens_maps");
        m_sm = ws.Get<CT>(m_smname);

		assert (params.exists("weights_name"));
		m_wname  = params.Get<std::string>("weights_name");
        
		if (params.exists("phase_cor"))
			m_pc = ws.Get<CT>(params.Get<std::string>("phase_cor"));
		if (params.exists("b0"))
			m_pc = ws.Get<T>(params.Get<std::string>("b0"));

		m_nx.push_back(ndims(m_sm)-1);
		Matrix<size_t> ms (m_nx[0],1);
		for (size_t i = 0; i < m_nx[0]; i++)
			ms[i] = size(m_sm,i);

        container<size_t> sizesm = vsize(m_sm);
        m_nx.push_back(sizesm.back()); // NC
        m_nx.push_back(numel(ws.Get<T>(m_wname))); // NK
        m_nx.push_back(std::accumulate(sizesm.begin(), sizesm.end(), 1, multiply<size_t>)/m_nx[1]); //NR

		m_cgiter  = params.Get<size_t>("cgiter");
		m_cgeps   = params.Get<double>("cgeps");
		m_lambda  = params.Get<double>("lambda");
		m_verbose = params.Get<int>("verbose");
		
#pragma omp parallel 
		{
		  m_np = omp_get_num_threads ();
		}
		
		if (params.exists("np") && boost::any_cast<int>(params["np"]) > 0) {
		  omp_set_num_threads(boost::any_cast<int>(params["np"]));
		}

		printf ("  Initialising NCSENSE:\n");
		printf ("  No of threads: %i\n", m_np);
		printf ("  Signal nodes: %li\n", m_nx[2]);
        printf ("  Channels: %zu\n", m_nx[1]);
        printf ("  Space size: %zu\n", m_nx[3]);
		printf ("  CG: eps(%.3e) iter(%li) lambda(%.3e)\n", m_cgeps, m_cgiter, m_lambda);
		printf ("  FT: eps(%.3e) iter(%li) m(%li) alpha(%.3e)\n", fteps, ftiter, m, alpha);

		for (size_t i = 0; i < m_np; i++) // FFTW planning not thread-safe
			m_fts.push_back(new NFFT<T> (ms, m_nx[2], m, alpha, b0, m_pc, fteps, ftiter));
		
		m_ic     = IntensityMap (m_sm);
		m_initialised = true;
        
		printf ("  ...done.\n\n");
		
	}

	/**
	 * @brief        Clean up and destruct NFFT plans
	 */ 
	virtual ~NCSENSE () {
		
		if (m_initialised)
			for (size_t i = 0; i < m_np; i++)
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
            m_fts[omp_get_thread_num()]->KSpace(k);
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
            m_fts[omp_get_thread_num()]->Weights(w);
		}
		
	}


	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<CT>
	Trafo       (const Matrix<CT>& m) const {

		return E (m, m_sm, m_nx, m_fts);

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
	virtual Matrix<CT>
	Trafo       (const Matrix<CT>& m, const Matrix<CT>& sens, const bool& recal = true) const {

		return E (m, sens, m_nx, m_fts);

	}

	
	/**
	 * @brief Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<CT>
	Adjoint     (const Matrix<CT>& m) const {

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
	virtual Matrix<CT>
	Adjoint (const Matrix<CT>& m,
			 const Matrix<CT>& sens,
			 const bool recal = false) const {

        T rn, rno, xn, ts;
		Matrix<CT> p, r, x, q;
		vector<T> res;
        std::vector< Matrix<cxfl> > vc;

		p  = EH (m, sens, m_nx, m_fts) * m_ic;
		r  = p;
		x  = zeros<CT>(size(p));
        xn = creal(p.dotc(p));
        rn = xn;
        if (m_verbose)
            vc.push_back (p/m_ic);

		for (size_t i = 0; i < m_cgiter; i++) {
			
			res.push_back(rn/xn);
			if (std::isnan(res.at(i)) || res.at(i) <= m_cgeps)  break;
 			printf ("    %03lu %.7f\n", i, res.at(i)); fflush (stdout);

			q  = EH(E(p * m_ic, sens, m_nx, m_fts), sens, m_nx, m_fts) * m_ic;
			if (m_lambda)
				q  += m_lambda * p;
            ts  = creal(rn);
            ts /= creal(p.dotc(q));
			x  += ts * p;
			r  -= ts * q;
			rno = rn;
			rn  = creal(r.dotc(r));
			p  *= rn / rno;
			p  += r;
            if (m_verbose)
                vc.push_back (x * m_ic);

		}

        if (m_verbose) {
            size_t cpsz = numel(x);
            x = Matrix<CT> (size(x,0), size(x,1), (m_nx[0] == 3) ? size(x,2) : 1, vc.size());
            typename container<CT>::iterator ti = x.Begin();
            
            for (size_t i = 0; i < vc.size(); i++) {
                std::copy (vc[i].Begin(), vc[i].End(), ti);
                ti += cpsz;
            }
        } else
            x *= m_ic;

        return x;

	}
	
	
	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<CT>
	operator* (const Matrix<CT>& m) const {
		return Trafo(m);
	}
	

	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<CT>
	operator->* (const Matrix<CT>& m) const {
		return Adjoint (m);
	}

	
	
private:

    std::vector<NFFT<T>*> m_fts;         /**< Non-Cartesian FT operators (Multi-Core?) */
	bool      m_initialised; /**< All initialised? */
    bool      m_verbose;

	Matrix<CT> m_sm;          /**< Sensitivities */
	Matrix<T> m_ic;     /**< Intensity correction I(r) */
	Matrix<CT> m_pc; /**< @brief Correction phase */

	std::string m_smname;
	std::string m_wname;

    std::vector<size_t> m_nx;
    
	size_t    m_cgiter;         /**< Max # CG iterations */
	double    m_cgeps;          /**< Convergence limit */
	double    m_lambda;         /**< Tikhonov weight */
	
	int       m_np;

};


#endif

