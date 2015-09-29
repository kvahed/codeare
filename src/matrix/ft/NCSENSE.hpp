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

enum FT_EXCEPTION {NCSENSE_KSPACE_DIMENSIONS};

template <class T> class NCSENSE : public FT<T>{

	// TODO: Check if k-space and weights have been assigned

    typedef typename TypeTraits<T>::RT RT;
    typedef Range<false> R;
    typedef Range<true> CR;
	
public:

	/**
	 * @brief         Default constructor
	 */
	NCSENSE() NOEXCEPT : m_initialised (false), m_cgiter(30), m_cgeps (1.0e-6), m_lambda (1.0e-6),
        m_verbose (false), m_np(0), m_3rd_dim_cart(false), m_nmany(1), m_dim4(1), m_dim5(1) {}
    
    
	/**
	 * @brief          Construct with parameters
	 *
	 * @param  params  Configuration parameters
	 */
	NCSENSE        (const Params& params) NOEXCEPT
              : FT<T>::FT(params), m_cgiter(30), m_initialised(false), m_cgeps(1.0e-6),
                m_lambda(1.0e-6), m_verbose (false), m_np(0), m_3rd_dim_cart(false), m_nmany(1), m_dim4(1), m_dim5(1)  {

		size_t cart_dim = 1;

		ft_params["epsilon"] = (params.exists("fteps")) ? fp_cast(params["fteps"]): (RT)1.0e-3;
		ft_params["alpha"] = params.exists("alpha") ? fp_cast(params["alpha"]) : (RT)1.;
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
		m_nx.push_back(unsigned_cast(params["nk"])); // NK
        m_nx.push_back(std::accumulate(sizesm.begin(), sizesm.end(), 1,
        		c_multiply<size_t>) / m_nx[1]); //NR

		m_cgiter  = params.Get<size_t>("cgiter");
		m_cgeps   = params.Get<RT>("cgeps");
		m_lambda  = params.Get<RT>("lambda");
        try {
            m_np  = params.Get<int>("threads");
        } catch (const boost::bad_any_cast&) {
            m_np  = std::thread::hardware_concurrency();
        }

        omp_set_num_threads(m_np);
        ft_params["threads"] = m_np;
        
        try {
            m_dim4  = params.Get<int>("dim4");
            m_dim5  = params.Get<int>("dim5");
            m_nmany = m_dim4*m_dim5;
        } catch (PARAMETER_MAP_EXCEPTION) {
        } catch (const boost::bad_any_cast&) {}
        
        try {
        	m_verbose = (params.Get<int>("verbose") > 0);
        } catch (const boost::bad_any_cast&) {}
        m_nx.push_back(m_np);
        
		ft_params["imsz"] = ms;

        if (m_nmany > 1) {
            for (size_t i = 0; i < m_nmany; ++i)
                m_fts.push_back(NFFT<T>(ft_params));
        } else{
            for (size_t i = 0; i < m_nx[1]; ++i)
                m_fts.push_back(NFFT<T>(ft_params));
        }
        
        omp_set_num_threads(m_fts.size());
        
		m_ic     = IntensityMap (m_sm);
		m_initialised = true;

        m_fwd_out = Matrix<T> (m_nx[2],cart_dim,m_nx[1],m_dim4,m_dim5); // nodes x slices x channels
        Vector<size_t> tmp = size(m_sm);
        if (m_nmany > 1) {
            tmp.PushBack(m_dim4);
            tmp.PushBack(m_dim5);
        }
		m_bwd_out = Matrix<T> (tmp);               // size of sensitivity maps
		
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
	void KSpace (const Matrix<RT>& k) {
		m_k = k;
        if (size(k,1) == KSpaceSize() && m_nmany == 1) {
#pragma omp parallel num_threads (m_fts.size())
            {
                m_fts[omp_get_thread_num()].KSpace(k);
            }
        } else if (size(m_k,2)*size(m_k,3) == m_nmany) {
#pragma omp parallel num_threads (m_fts.size())
            {
                size_t i = omp_get_thread_num(), l=i%size(m_k,2), n = i/size(m_k,2);
                if (ndims(k)==4)
                    m_fts[i].KSpace(k(CR(),CR(),CR(l),CR(n)));
                else if (ndims(k) == 5)
                    m_fts[i].KSpace(k(CR(),CR(),CR(),CR(l),CR(n)));
                else 
                    throw NCSENSE_KSPACE_DIMENSIONS;
            }
        } else {
            throw NCSENSE_KSPACE_DIMENSIONS;
        }
                                
    }
	
    
	/**
	 * @brief      Assign k-space weigths (jacobian of k in t) 
	 * 
	 * @param  w   Weights
	 */
	void Weights (const Matrix<RT>& w) NOEXCEPT {
		m_w = w;
        for (size_t i = 0; i < m_fts.size(); ++i)
            m_fts[i].Weights(w);
	}
    
    
	virtual Matrix<T> operator/ (const MatrixType<T>& m) const NOEXCEPT {
        Matrix<T> ret;
        if (m_nmany > 1) {
#pragma omp parallel num_threads (m_nmany)
            {
                size_t k = omp_get_thread_num(), l = k%m_dim4, n = k/m_dim4;
                for (int j = 0; j < m_nx[1]; ++j)
                    if (m_nx[0] == 2)
                        m_bwd_out (R(),R(),    R(j),R(l),R(n)) =
                            m_fts[k] ->* m(CR(),     CR(j),CR(l),CR(n));
                    else
                        m_bwd_out (R(),R(),R(),R(j),R(l),R(n)) =
                            m_fts[k] ->* m(CR(),CR(),CR(j),CR(l),CR(n));
                m_bwd_out(R(),R(),R(),R(),R(l),R(n)) *= m_csm;
            }
            ret = squeeze(sum(m_bwd_out,3));
        } else {
#pragma omp parallel num_threads (m_nx[1])
            {
                size_t k = omp_get_thread_num();
                if (m_nx[0] == 2)
                    m_bwd_out (R(),R(),    R(k)) = m_fts[k] ->* m(CR(),     CR(k));
                else
                    m_bwd_out (R(),R(),R(),R(k)) = m_fts[k] ->* m(CR(),CR(),CR(k));
            }
            ret = squeeze(sum(m_bwd_out*m_csm,size(m_sm).size()-1)) * m_ic;
        }
	    return ret;
	}
    
    
	/**
  	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<T> Trafo (const MatrixType<T>& m) const NOEXCEPT {
        if (m_nmany > 1) {
#pragma omp parallel num_threads (m_nmany)
            {
                size_t k = omp_get_thread_num(), l = k%m_dim4, n = k/m_dim4;
                if (m_nx[0] == 2)
                    for (int j = 0; j < m_nx[1]; ++j)
                        m_fwd_out(R(),    R(j),R(l),R(n)) =
                            m_fts[k] * (m_sm(CR(),CR(),     CR(j))*m(CR(),CR(),     CR(l),CR(n)));
                else
                    for (int j = 0; j < m_nx[1]; ++j)
                        m_fwd_out(R(),R(),R(j),R(l),R(n)) =
                            m_fts[k] * (m_sm(CR(),CR(),CR(),CR(j))*m(CR(),CR(),CR(),CR(l),CR(n)));
            }            
        } else {
#pragma omp parallel num_threads (m_fts.size())
            {
                size_t j = omp_get_thread_num();
                if (m_3rd_dim_cart)
                    m_fwd_out(R(),R(),R(),R(j)) = m_fts[j] * (m_sm(CR(),CR(),CR(),CR(j))*m);
                else
                    m_fwd_out(R(),R(),    R(j)) = m_fts[j] * (m_sm(CR(),CR(),     CR(j))*m);
            }
        }
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
	virtual Matrix<T> Trafo (const MatrixType<T>& m, const MatrixType<T>& sens,
                             const bool& recal = true) const NOEXCEPT {
		return Trafo(m);
	}
    
	
	/**
	 * @brief Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<T> Adjoint (const MatrixType<T>& m) const NOEXCEPT {
		return this->Adjoint (m, m_sm, false);
	}
	
	
	/**
	 * @brief Estimate coil sensitivities
	 * @param  data  measurement
	 */
	virtual void EstimateSensitivities (const MatrixType<T>& data, size_t nk) const {
		Matrix<T> out (size(m_sm));
#pragma omp parallel for
		for (int i = 0; i < m_nx[1]; ++i)
			m_fts[omp_get_thread_num()].NFFTPlan().M_total = nk;
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
	virtual Matrix<T> Adjoint (const MatrixType<T>& m, const MatrixType<T>& sens,
                               const bool recal = false) const NOEXCEPT {
		// TODO: Not functional yet
		if (m_sm.Size() == 1)
			EstimateSensitivities(m, m_nx[2]);
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
	
    
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	virtual Matrix<T> operator->* (const MatrixType<T>& m) const NOEXCEPT {
		return Adjoint (m);
	}
    
	virtual std::ostream& Print (std::ostream& os) const {
		Operator<T>::Print(os);
		os << "    NCCG: eps("<< m_cgeps << ") iter(" << m_cgiter << ") lambda(" << m_lambda << ")" << std::endl;
		os << "    threads(" << m_np << ") channels(" << m_nx[1] << ") nmany(" << m_nmany << ")" << std::endl;
		os << m_fts[0];
		return os;
	}

    virtual FT<T>* getFT () {
        return (FT<T>*) &m_fts[0];
    }

    inline size_t KSpaceSize () const {
        return m_fts[0].KSpaceSize();
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
	RT     m_cgeps;          /**< Convergence limit */
	RT     m_lambda;         /**< Tikhonov weight */

    size_t     m_nmany;          /**< Recounstruct multiple volumes with same k-space and maps */
	size_t m_dim4, m_dim5;
	int        m_np;

    

	Params ft_params;
	Matrix<RT> m_k, m_w;

	mutable Matrix<T> m_fwd_out, m_bwd_out;

	mutable codeare::optimisation::CGLS<T> m_cgls;

};


#endif

