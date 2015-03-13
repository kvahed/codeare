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
#include "mri/MRI.hpp"
#include "Lapack.hpp"
#include "tinyxml.h"

#include "Workspace.hpp"

#include <numeric>
#include <boost/math/special_functions/fpclassify.hpp>

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
	NCSENSE() NOEXCEPT : m_initialised (false),
                m_cgiter(30),
                m_cgeps (1.0e-6),
                m_lambda (1.0e-6),
                m_verbose (false),
                m_np(0),
                m_3rd_dim_cart(false) {}
    
    
	/**
	 * @brief          Construct with parameters
	 *
	 * @param  params  Configuration parameters
	 */
	NCSENSE        (const Params& params) NOEXCEPT 
              : FT<T>::FT(params),
                m_cgiter(30),
                m_initialised(false),
                m_cgeps(1.0e-6),
                m_lambda(1.0e-6),
                m_verbose (false),
                m_np(0),
                m_3rd_dim_cart(false) {

		Params ft_params;
		size_t cart_dim = 1;

		if (params.exists("fteps"))
			ft_params["epsilon"] = fp_cast(params["fteps"]);
		if (params.exists("alpha"))
			ft_params["alpha"] = fp_cast(params["alpha"]);
		if (params.exists("ftiter"))
			ft_params["maxit"] = unsigned_cast(params["ftiter"]);
		if (params.exists("m"))
			ft_params["m"] = unsigned_cast(params["m"]);
		if (params.exists("nk"))
			ft_params["nk"] = unsigned_cast(params["nk"]);

		if (params.exists("3rd_dim_cart")) {
			try {
				m_3rd_dim_cart = params.Get<bool>("3rd_dim_cart");
			} catch (const boost::bad_any_cast&) {
				printf ("  WARNING - NCSENSE: Could not interpret input for Cartesian "
						"nature of 3rd dimension. \n");
			}
		}
		ft_params["3rd_dim_cart"] = m_3rd_dim_cart;

		Workspace& ws = Workspace::Instance();
		Matrix<T> b0;

		assert (params.exists("sensitivities"));
		m_sm = params.Get<Matrix<CT> >("sensitivities");
		if (m_3rd_dim_cart)
			cart_dim = size(m_sm,2);

		if (params.exists("phase_cor"))
			m_pc = ws.Get<CT>(params.Get<std::string>("phase_cor"));
		if (params.exists("b0"))
			m_pc = ws.Get<T>(params.Get<std::string>("b0"));
		
		m_nx.push_back(ndims(m_sm)-1);
		Vector<size_t> ms (m_nx[0]);
		for (size_t i = 0; i < m_nx[0]; i++)
			ms[i] = size(m_sm,i);

		Vector<size_t> sizesm = size(m_sm);
        m_nx.push_back(sizesm.back());             // NC
		m_nx.push_back(unsigned_cast(params["nk"])*cart_dim); // NK
        m_nx.push_back(std::accumulate(sizesm.begin(), sizesm.end(), 1, c_multiply<size_t>)/m_nx[1]); //NR

		m_cgiter  = params.Get<size_t>("cgiter");
		m_cgeps   = params.Get<double>("cgeps");
		m_lambda  = params.Get<double>("lambda");
		m_verbose = (params.Get<int>("verbose") > 0);
			
		printf ("  Initialising NCSENSE:\n");
		printf ("  No of threads: %i\n", m_np); // TODO: Thread num not passed on
		printf ("  Signal nodes: %li\n", m_nx[2]);
		printf ("  Cartesian 3rd dimension: %d\n", m_3rd_dim_cart);
        printf ("  Channels: " JL_SIZE_T_SPECIFIER "\n", m_nx[1]);
        printf ("  Space size: " JL_SIZE_T_SPECIFIER "\n", m_nx[3]);
		printf ("  CG: eps(%.3e) iter(%li) lambda(%.3e)\n", m_cgeps, m_cgiter, m_lambda);
		//printf ("  FT: eps(%.3e) iter(%li) m(%li) alpha(%.3e)\n", fteps, ftiter, m, alpha);

		ft_params["imsz"] = ms;

        m_fts = NFFT<T>(ft_params);

		m_ic     = IntensityMap (m_sm);
		m_initialised = true;

		printf ("  ...done.\n\n");
		
	}


	/**
	 * @brief        Clean up and destruct NFFT plans
	 */ 
	virtual ~NCSENSE () NOEXCEPT {}
	

	/**
	 * @brief      Assign k-space trajectory
	 * 
	 * @param  k   K-space trajectory
	 */
	void
	KSpace (const Matrix<T>& k) NOEXCEPT {
        m_fts.KSpace(k);
	}
	

	/**
	 * @brief      Assign k-space weigths (jacobian of k in t) 
	 * 
	 * @param  w   Weights
	 */
	void
	Weights (const Matrix<T>& w) NOEXCEPT {		
        m_fts.Weights(w);
	}


	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<CT>
	Trafo       (const Matrix<CT>& m) const NOEXCEPT {
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
	Trafo       (const Matrix<CT>& m, const Matrix<CT>& sens, const bool& recal = true) const NOEXCEPT {
		return E (m, sens, m_nx, m_fts);
	}

	
	/**
	 * @brief Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<CT>
	Adjoint     (const Matrix<CT>& m) const NOEXCEPT {
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
			 const bool recal = false) const NOEXCEPT {

        T rn, rno, xn, ts;
		Matrix<CT> p, r, x, q;
		vector<T> res;
        Vector<Matrix<cxfl> > vc;

        typedef typename Vector<CT>::iterator it_type;

		p = EH (m, sens, m_nx, m_fts)* m_ic;
		if (m_cgiter == 0)
			return p;

		r  = p;
        xn = real(p.dotc(p));
        rn = xn;

        if (m_verbose)
            vc.push_back (p);

		for (size_t i = 0; i < m_cgiter; i++) {
			res.push_back(rn/xn);
			if (i==0)
				x  = zeros<CT>(size(p));
			if (boost::math::isnan(res.at(i)) || res.at(i) <= m_cgeps)
				break;
			if (m_verbose)
				printf ("    %03lu %.7f\n", i, res.at(i));
			q  = EH(E(p * m_ic, sens, m_nx, m_fts), sens, m_nx, m_fts) * m_ic;

			if (m_lambda)
				q  += m_lambda * p;
			ts  = rn / real(p.dotc(q));
			x  += ts * p;
			r  -= ts * q;
			rno = rn;
			rn  = real(r.dotc(r));
			p  *= rn / rno;
			p  += r;
			if (m_verbose)
				vc.push_back(x * m_ic);
		}

        if (m_verbose) { // Keep intermediate results
            size_t cpsz = numel(x);
            x = Matrix<CT> (size(x,0), size(x,1), (m_nx[0] == 3) ? size(x,2) : 1, vc.size());
            it_type it = x.Begin();
            for (size_t i = 0; i < vc.size(); i++) {
                std::copy (vc[i].Begin(), vc[i].End(), it);
                it += cpsz;
            }
            vc.Clear();
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
	operator* (const Matrix<CT>& m) const NOEXCEPT {
		return Trafo(m);
	}
	

	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<CT>
	operator->* (const Matrix<CT>& m) const NOEXCEPT {
		return Adjoint (m);
	}

	
	
private:

    //Vector<NFFT<T> > m_fts; /**< Non-Cartesian FT operators (Multi-Core?) */
    NFFT<T>    m_fts;
	bool       m_initialised; /**< All initialised? */
    bool       m_verbose;	  /**< Verbose binary output (keep all intermediate steps) */
    bool       m_3rd_dim_cart; /**< 3rd FT dimension is Cartesian (stack of ...) */

	Matrix<CT> m_sm;          /**< Sensitivities */
	Matrix<T>  m_ic;     /**< Intensity correction I(r) */
	Matrix<CT> m_pc; /**< @brief Correction phase */

	std::string m_smname;     /**< Sensitivity map name */
	std::string m_wname;	  /**< Weights name */

    Vector<size_t> m_nx;
    
	size_t     m_cgiter;         /**< Max # CG iterations */
	double     m_cgeps;          /**< Convergence limit */
	double     m_lambda;         /**< Tikhonov weight */
	
	int        m_np;

};


#endif

