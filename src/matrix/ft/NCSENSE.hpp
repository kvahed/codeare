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
//#include "SEM.hpp"
#include "mri/MRI.hpp"
#include "Lapack.hpp"
#include "tinyxml.h"
#include "IOContext.hpp"
#include "CGLS.hpp"

#include "Workspace.hpp"

#include <numeric>
#include <thread>

#include <boost/math/special_functions/fpclassify.hpp>

/**
 * @brief Non-Cartesian SENSE<br/>
 *        According Pruessmann et al. (2001). MRM, 46(4), 638-51.
 *
 */
template <class T>
class NCSENSE : public FT<T>{

	// TODO: Check if k-space and weights have been assigned

    typedef typename TypeTraits<T>::RT RT;
    typedef Range<false> R;
    typedef Range<true> CR;
	
public:

	/**
	 * @brief         Default constructor
	 */
	NCSENSE() NOEXCEPT : m_initialised (false), m_cgiter(30), m_cgeps (1.0e-6), m_lambda (1.0e-6),
                m_verbose (false), m_np(0), m_3rd_dim_cart(false) {}
    
    
	/**
	 * @brief          Construct with parameters
	 *
	 * @param  params  Configuration parameters
	 */
	NCSENSE        (const Params& params) NOEXCEPT
              : FT<T>::FT(params), m_cgiter(30), m_initialised(false), m_cgeps(1.0e-6),
                m_lambda(1.0e-6), m_verbose (false), m_np(0), m_3rd_dim_cart(false) {

		size_t cart_dim = 1;

		ft_params["epsilon"] = (params.exists("fteps")) ? fp_cast(params["fteps"]): 1.0e-3;
		ft_params["alpha"] = params.exists("alpha") ? fp_cast(params["alpha"]) : 1.;
		ft_params["maxit"] = (params.exists("ftiter")) ? unsigned_cast(params["ftiter"]) : 3;
		ft_params["m"] = (params.exists("m")) ? unsigned_cast(params["m"]) : 1;
		ft_params["nk"] = (params.exists("nk")) ? unsigned_cast(params["nk"]) : 1;

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
		Matrix<RT> b0;

		assert (params.exists("sensitivities"));
		m_sm = params.Get<Matrix<T> >("sensitivities");
        m_csm = conj(m_sm);
		if (m_3rd_dim_cart)
			cart_dim = size(m_sm,2);

		if (params.exists("phase_cor"))
			m_pc = ws.Get<T>(params.Get<std::string>("phase_cor"));
		if (params.exists("b0"))
			m_pc = ws.Get<RT>(params.Get<std::string>("b0"));
		
		m_nx.push_back(ndims(m_sm)-1);
		Vector<size_t> ms (m_nx[0]);
		for (size_t i = 0; i < m_nx[0]; i++)
			ms[i] = size(m_sm,i);

		Vector<size_t> sizesm = size(m_sm);
        m_nx.push_back(sizesm.back());             // NC
		m_nx.push_back(unsigned_cast(params["nk"])*cart_dim); // NK
        m_nx.push_back(std::accumulate(sizesm.begin(), sizesm.end(), 1,
        		c_multiply<size_t>) / m_nx[1]); //NR

		m_cgiter  = params.Get<size_t>("cgiter");
		m_cgeps   = params.Get<double>("cgeps");
		m_lambda  = params.Get<double>("lambda");
        try {
            m_np  = params.Get<int>("threads");
        } catch (const boost::bad_any_cast&) {
            m_np  = std::thread::hardware_concurrency();
        }
        omp_set_num_threads(m_np);
        ft_params["threads"] = m_np;

        try {
        	m_verbose = (params.Get<int>("verbose") > 0);
        } catch (const boost::bad_any_cast&) {}
        m_nx.push_back(m_np);
        
		ft_params["imsz"] = ms;

        for (size_t i = 0; i < m_np; ++i)
            m_fts.PushBack(NFFT<T>(ft_params));

		m_ic     = IntensityMap (m_sm);
		m_initialised = true;

		if (m_3rd_dim_cart)
            m_fwd_out = Matrix<T> (m_nx[2]/cart_dim,cart_dim,m_nx[1]);
        else
            m_fwd_out = Matrix<T> (m_nx[2],cart_dim,m_nx[1]);
		m_bwd_out = Matrix<T> (size(m_sm));
		
		m_cgls = codeare::optimisation::CGLS<T>(m_cgiter, m_cgeps, m_lambda, m_verbose);

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
	KSpace (const Matrix<RT>& k) NOEXCEPT {
		m_k = k;
        for (size_t i = 0; i < m_fts.size(); ++i)
            m_fts[i].KSpace(k);
	}
	

	/**
	 * @brief      Assign k-space weigths (jacobian of k in t) 
	 * 
	 * @param  w   Weights
	 */
	void
	Weights (const Matrix<RT>& w) NOEXCEPT {
		m_w = w;
        for (size_t i = 0; i < m_fts.size(); ++i)
            m_fts[i].Weights(w);
	}


	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<T>
	Trafo       (const MatrixType<T>& m) const NOEXCEPT {
#pragma omp parallel for
	    for (int j = 0; j < m_nx[1]; ++j)
            if (m_3rd_dim_cart)
                0;//Slice  (m_fwd_out, j, m_fts[k] * (resize(Volume (m_sm, j),size(m)) * m));
            else
                m_fwd_out(R(),R(),R(j)) = m_fts[omp_get_thread_num()] * (m_sm(CR(),CR(),CR(j))*m);
//                    (m /* * ((m_nx[0] == 2) ? m_sm(CR(),CR(),CR(j)) : m_sm(CR(),CR(),CR(),CR(j)))*/);
	    return squeeze(m_fwd_out);
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
	virtual Matrix<T>
	Trafo       (const MatrixType<T>& m, const MatrixType<T>& sens,
			const bool& recal = true) const NOEXCEPT {
		return Trafo(m);
	}

	
	/**
	 * @brief Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<T>
	Adjoint     (const MatrixType<T>& m) const NOEXCEPT {
		return this->Adjoint (m, m_sm, false);
	}
	
	
	/**
	 * @brief Estimate coil sensitivities
	 * @param  data  measurement
	 */
	virtual void
	EstimateSensitivities (const MatrixType<T>& data, size_t nk = 8192) const {
		Matrix<T> out (size(m_sm));
#pragma omp parallel for
		for (int i = 0; i < m_nx[1]; ++i) {
			m_fts[omp_get_thread_num()].NFFTPlan().M_total = nk;
/*			if (m_nx[0] == 2)
				Slice  (out, i,
						m_fts[omp_get_thread_num()] ->* data(CR(0,8191), CR(i)));
			else
				Volume (out, i,
                m_fts[omp_get_thread_num()] ->* data(CR(0,8191), CR(i)));*/
            m_fts[omp_get_thread_num()].NFFTPlan().M_total = m_nx[2];
		}
	}


	/**
	 * @brief Image size
	 * @return image size
	 */
	inline Vector<size_t> ImageSize () {
		Vector<size_t> ret;
		//TODO: Image size, K-space size
		return ret;
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
	virtual Matrix<T>
	Adjoint (const MatrixType<T>& m,
			 const MatrixType<T>& sens,
			 const bool recal = false) const NOEXCEPT {

		// TODO: Not functional yet
		if (m_sm.Size() == 1)
			EstimateSensitivities(m);
        return m_cgls.Solve(*this, m);

	}
	
	
	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<T> operator* (const MatrixType<T>& m) const NOEXCEPT {
		return Trafo(m);
	}
	

	virtual Matrix<T> operator/ (const MatrixType<T>& m) const NOEXCEPT {
#pragma omp parallel for
		for (int j = 0; j < m_nx[1]; ++j) 
	        if (m_nx[0] == 2)
	        	m_bwd_out (R(),R(),    R(j)) =
                    m_fts[omp_get_thread_num()] ->* m(CR(),     CR(j));
	        else
	        	m_bwd_out (R(),R(),R(),R(j)) =
                    m_fts[omp_get_thread_num()] ->* m(CR(),CR(),CR(j));
	    return squeeze(sum(m_bwd_out*m_csm,size(m_sm).size()-1)) * m_ic;
	}


	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<T>
	operator->* (const MatrixType<T>& m) const NOEXCEPT {
		return Adjoint (m);
	}

	virtual std::ostream& Print (std::ostream& os) const {
		FT<T>::Print(os);
		os << "    NCCG: eps("<< m_cgeps << ") iter(" << m_cgiter <<
				") lambda(" << m_lambda << ")" << std::endl;
		os << "    threads(" << m_np << ") channels(" << m_nx[1] << ")" << std::endl;
		os << m_fts[0];
		return os;
	}
	
private:

	mutable Vector<NFFT<T> > m_fts; /**< Non-Cartesian FT operators (Multi-Core?) */
	bool       m_initialised; /**< All initialised? */
    bool       m_verbose;	  /**< Verbose binary output (keep all intermediate steps) */
    bool       m_3rd_dim_cart; /**< 3rd FT dimension is Cartesian (stack of ...) */

	mutable Matrix<T> m_sm;          /**< Sensitivities */
	mutable Matrix<T> m_csm;          /**< Sensitivities */
	Matrix<RT>  m_ic;     /**< Intensity correction I(r) */
	Matrix<T> m_pc; /**< @brief Correction phase */

	std::string m_smname;     /**< Sensitivity map name */
	std::string m_wname;	  /**< Weights name */

    Vector<size_t> m_nx;
    
	size_t     m_cgiter;         /**< Max # CG iterations */
	double     m_cgeps;          /**< Convergence limit */
	double     m_lambda;         /**< Tikhonov weight */
	
	int        m_np;

	Params ft_params;
	Matrix<RT> m_k, m_w;

	mutable Matrix<T> m_fwd_out, m_bwd_out;

	mutable codeare::optimisation::CGLS<T> m_cgls;

};


#endif

