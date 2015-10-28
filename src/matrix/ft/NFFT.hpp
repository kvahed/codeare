/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
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

#include "NFFTTraits.hpp"
#include "FFTWTraits.hpp"
#include "Algos.hpp"
#include "Access.hpp"
#include "FT.hpp"
#include "CX.hpp"
#include "Creators.hpp"

#include <thread>

/**
 * @brief Matrix templated ND non-equidistand Fourier transform with NFFT 3 (TU Chemnitz)<br/>
 */
template <class T>
class NFFT : public FT<T> {

#ifdef HAVE_NFFT3F_UNUSED
    typedef T NFFTType;
    typedef typename TypeTraits<T>::RT NFFTRType;
#else
    typedef std::complex<double> NFFTType;
    typedef double NFFTRType;
#endif
	typedef typename TypeTraits<T>::RT RT;
    typedef typename NFFTTraits<NFFTType>::Plan   Plan;
    typedef typename NFFTTraits<NFFTType>::B0Plan B0Plan;
    typedef typename NFFTTraits<NFFTType>::Solver Solver;
    typedef typename FTTraits<T>::Plan CartPlan;

public:

    /**
     * @brief         Default constructor
     */
    NFFT() NOEXCEPT :  m_initialised (false), m_have_pc (false), m_imgsz (0),
        m_M (0), m_maxit (0), m_rank (0), m_m(0), m_have_b0(false),
        m_3rd_dim_cart(false),m_ncart(1), m_alpha(1.), m_have_weights(false),
		m_have_kspace(false), m_np(std::thread::hardware_concurrency()),
        m_per_slice_kspace(false) {};

    /**
     * @brief        Construct with parameter set
     */
    inline NFFT (const Params& p) NOEXCEPT : m_have_b0(false), m_3rd_dim_cart(false),
        m_t (Matrix<RT>()), m_b0 (Matrix<RT>()), m_maxit(3), m_m(1), m_alpha(1.),
		m_epsilon(7.e-4f), m_sigma(1.0), m_ncart(1),  m_have_weights(false),
        m_have_kspace(false), m_np(std::thread::hardware_concurrency()),
        m_per_slice_kspace(false) {

        if (p.exists("nk")) {// Number of kspace samples
            try {
                m_M = unsigned_cast(p["nk"]);
            } catch (const boost::bad_any_cast& e) {
                printf ("**ERROR - NFFT: Numer of ksppace samples need to be "
                		"specified\n%s\n", e.what());
                assert(false);
            }
        } else {
            printf ("**ERROR - NFFT: Numer of ksppace samples need to be specified\n");
            assert(false);
        }

        m_3rd_dim_cart = try_to_fetch (p, "3rd_dim_cart", false);

        if (p.exists("imsz")) {// Image domain size
            try {
                m_N = boost::any_cast<Vector<size_t> >(p["imsz"]);
            } catch (const boost::bad_any_cast& e) {
                printf ("**ERROR - NFFT: Image domain dimensions need to be "
                		"specified\n%s\n", e.what());
                assert(false);
            }
        } else {
            printf ("**ERROR - NFFT: Image domain dimensions need to be specified\n");
            assert(false);
        }

        if (m_3rd_dim_cart) { // 3rd dimension is Cartesian
        	m_ncart = m_N.back();
        	m_N.pop_back();
        }
        m_n       = m_N;
        
        m_m       = try_to_fetch<size_t>(p, "m", 1);
        m_alpha   = try_to_fetch(p, "alpha", 1.0f);
        for (size_t i = 0; i < m_N.size(); ++i)
            m_n[i] = ceil(m_alpha*m_N[i]);
        
        m_rank    = m_N.size();
        m_imgsz   = 2*prod(m_N);

        m_epsilon = try_to_fetch(p, "epsilon", 7.0e-4f);
        m_maxit   = try_to_fetch<size_t>(p, "maxit", 3);

        if (p.exists("b0")) {
            try {
                m_b0 = p.Get<Matrix<RT> >("b0");
                if (numel(m_b0) != prod(m_N))
                    printf ("  WARNING - NFFT: b0 field does not fit to image space. "
                            "Defaulting to flat b0 = 0\n");
                else
                    m_have_b0 = true;
            } catch (const boost::bad_any_cast&) {
                printf ("  WARNING - NFFT: Could not interpret input for b0. "
                        "Defaulting to flat b0 = 0\n");
            }
        }
        
        if (p.exists("timing")) {
            try {
                m_t = p.Get<Matrix<RT> >("timing");
                if (numel(m_t) != m_M)
                    printf ("  WARNING - NFFT: kspace trajectory timing fit "
                    		"to trajectory. Defaulting to flat b0 = 0\n");
                else
                    m_have_b0 = (m_have_b0 && true);
            } catch (const boost::bad_any_cast&) {
                printf ("  WARNING - NFFT: Could not interpret input for kspace "
                		"trajectory timing. Defaulting to flat b0 = 0\n");
            }
        }
        

        if (m_have_b0) { // b0
            
            m_min_t  = min(m_t.Container());
            m_max_t  = max(m_t.Container());
            m_min_b0 = min(m_b0.Container());
            m_max_b0 = max(m_b0.Container());
            m_sigma  = 1.2;
            
            m_N.push_back(std::ceil(std::max(
            		fabs(m_min_b0),fabs(m_max_b0)) * m_max_t-m_min_t/2. +
            		(m_m)/(2.*m_sigma))*4.*m_sigma);
            
            if (m_N.back()%2!=0)
                m_N.back()++; // need even dimension

            m_n.push_back(m_N.back());
            m_n *= m_sigma;

            m_w = std::max(std::abs(m_min_b0),std::abs(m_max_b0))/(.5-((RT) m_m)/m_N[2]);
            m_ts =  (m_min_t+m_max_t)/2.;
            RT t    = ((m_max_t-m_min_t)/2.)/(.5-((RT) (m_m))/m_N[2]);

        	NFFTTraits<NFFTType>::Init (m_N, m_M, m_n, m_m, m_sigma,
        			m_b0_plan, m_solver);

            for (size_t j = 0; j < m_N[0]*m_N[1]; ++j)
                m_b0_plan.w[j] = m_b0[j] / m_w;


        } else {
        	NFFTTraits<NFFTType>::Init (m_N, m_M, m_n, m_m, m_plan, m_solver);
        }
        
        if (p.exists("pc")) {
            try {
                m_pc = boost::any_cast<T>(p["pc"]);
            } catch (const boost::bad_any_cast& e) {
                printf ("  WARNING - NFFT: Phase correction matrix corrupt?\n%s\n", e.what());
            }
            m_cpc = conj(m_pc);
            m_have_pc = true;
        }

        m_initialised = true;
        
    }
    
    /**
     * @brief Copy conctructor
     */
    NFFT (const NFFT<T>& ft) NOEXCEPT {
        *this = ft;
    }

    
    /**
     * @brief        Clean up and destruct
     */ 
    virtual ~NFFT () NOEXCEPT {
        if (m_initialised)
        	if (m_have_b0)
        		NFFTTraits<NFFTType>::Finalize (m_b0_plan, m_solver);
        	else
        		NFFTTraits<NFFTType>::Finalize (m_plan, m_solver);
    }
    
    
    /**
     * @brief     Assignement
     */
    inline NFFT<T>& operator= (const NFFT<T>& ft) NOEXCEPT {
        m_initialised = ft.m_initialised;
        m_have_pc     = ft.m_have_pc;
        m_rank        = ft.m_rank;
        m_pc          = ft.m_pc;
        m_cpc         = ft.m_cpc;
        m_N           = ft.m_N;
        m_n           = ft.m_n;
        m_M           = ft.m_M;
        m_maxit       = ft.m_maxit;
        m_epsilon     = ft.m_epsilon;
        m_imgsz       = ft.m_imgsz;
        m_m           = ft.m_m;
        m_have_b0     = ft.m_have_b0;
        m_b0          = ft.m_b0;
        m_t           = ft.m_t;
        m_min_t       = ft.m_min_t;
        m_max_t       = ft.m_max_t;
        m_min_b0      = ft.m_min_b0;
        m_max_b0      = ft.m_max_b0;
        m_sigma       = ft.m_sigma;
        m_alpha       = ft.m_alpha;
        m_3rd_dim_cart = ft.m_3rd_dim_cart;
        m_ncart       = ft.m_ncart;
        m_np          = ft.m_np;
        m_per_slice_kspace = ft.m_per_slice_kspace;
        if (m_have_b0)
        	NFFTTraits<NFFTType>::Init (m_N, m_M, m_n, m_m, m_sigma, m_b0_plan, m_solver);
        else
        	NFFTTraits<NFFTType>::Init (m_N, m_M, m_n, m_m,          m_plan,    m_solver);
        return *this;
    }
    
    /**
     * @brief      Assign k-space 
     * 
     * @param  k   Kspace trajectory
     */
    inline virtual void KSpace (const Matrix<RT>& k) {
        m_k = k;
        if (m_have_b0) { // +1D for omega
            for (size_t j = 0; j < (size_t)m_b0_plan.M_total; ++j) {
				m_b0_plan.plan.x[3*j+0] = k[2*j+0];
				m_b0_plan.plan.x[3*j+1] = k[2*j+1];
                m_b0_plan.plan.x[3*j+2] = (m_t[j]-m_ts)*m_w/m_N.back();
            }
        } else {
            if (k.Size() == m_plan.M_total*m_rank)
                std::copy (k.Begin(), k.End(), m_plan.x);
            else if (k.Size() == m_plan.M_total*m_rank*m_ncart)
                m_per_slice_kspace = true;
        }
        m_have_kspace = true;
    }
    
    
    /**
     * @brief      Assign k-space weigths (jacobian of trajectory with regards to time) 
     * 
     * @param  w   Weights
     */
    inline virtual void Weights (const Matrix<RT>& w) NOEXCEPT {
        m_kw = w;
    	if (m_have_b0)
    		assert (w.Size() == m_b0_plan.M_total);
    	else
    		assert (w.Size() == m_plan.M_total);
        std::copy (w.Begin(), w.End(), m_solver.w);
        if (m_have_b0) {
            NFFTTraits<NFFTType>::Weights (m_b0_plan.plan, m_solver, m_rank);
            NFFTTraits<NFFTType>::Psi (m_b0_plan.plan);
        } else {
        	NFFTTraits<NFFTType>::Weights (m_plan, m_solver, m_rank);
        	NFFTTraits<NFFTType>::Psi (m_plan);
        }
        m_have_weights = true;
    }
    
    
    /**
     * @brief    Forward transform
     *
     * @param  m To transform
     * @return   Transform
     */
    inline virtual Matrix<T>
    Trafo       (const MatrixType<T>& m) const NOEXCEPT {

		NFFTRType* tmpd;
		RT* tmpt;
        Matrix<T> out (m_M, ((m_3rd_dim_cart && m_ncart > 1) ? m_ncart : 1));

/*        if (m_3rd_dim_cart && m_ncart > 1) { // Cartesian FT 3rd dim
        	int n = static_cast<int>(m_ncart);
        	tmpm = permute (tmpm, 2, 0, 1);
        	size_t cent = std::ceil(.5*m_ncart);
        	Matrix<T> tmp;
        	for (size_t i = 0; i < m_imgsz/2; ++i) { // fftshift colums
				tmp = Column (tmpm,i);
				std::rotate(tmp.Begin(), tmp.Begin()+cent, tmp.End());
				Column (tmpm, i, tmp);
        	}
        	CartPlan cp = FTTraits<T>::DFTPlanMany (1, &n, m_imgsz/2,
        			(FTType*)tmpm.Ptr(), (FTType*)tmpm.Ptr(), FFTW_FORWARD);
        	FTTraits<T>::Execute(cp);
        	FTTraits<T>::Destroy(cp);
        	for (size_t i = 0; i < m_imgsz/2; ++i) { // fftshift colums
				tmp = Column (tmpm,i);
				std::rotate(tmp.Begin(), tmp.Begin()+cent, tmp.End());
				Column (tmpm, i, tmp);
        	}
        	tmpm = permute (tmpm, 1, 2, 0)/sqrt((RT)m_ncart);
            }*/

        size_t tmp = numel(m)/m_ncart;
        for (size_t i = 0; i < m_ncart; ++i) {

			tmpd = (NFFTRType*) m_plan.f_hat;
            size_t os = i*tmp;
                
            for (size_t j = 0; j < tmp; ++j) {

                T val = (!m_have_b0) ? m[j+os] : m[j+os] *
                    std::polar<RT>((RT)1., (RT)(2. * PI * m_ts * m_b0[j] * m_w));
                tmpd[2*j+0] = real(m[j+os]);
                tmpd[2*j+1] = imag(m[j+os]);
            }

			if (m_have_b0)
				NFFTTraits<NFFTType>::Trafo (m_b0_plan);
			else
				NFFTTraits<NFFTType>::Trafo (m_plan);

			tmpd = (NFFTRType*) m_plan.f;
			tmpt = (RT*) (out.Ptr() + i*m_M);
			std::copy (tmpd, tmpd+2*m_M, tmpt);

        }

        return squeeze(out);
        
    }
    
    
    /**
     * @brief    Backward transform
     *
     * @param  m To transform
     * @return   Transform
     */
	virtual Matrix<T> Adjoint (const MatrixType<T>& m) const {

        Vector<size_t> N = m_N;
        NFFTRType* tmpd;
        RT* tmpt;
        
        if (m_have_b0)
            N.pop_back();
        if (m_3rd_dim_cart && m_ncart > 1) // Cartesian FT 3rd dim
        	N.push_back(m_ncart);

        Matrix<T> out (N);
        for (size_t i = 0; i < m_ncart; ++i) {
            
			tmpd = (NFFTRType*) m_solver.y;
            
            size_t os = i*m_M;
            for (size_t j = 0; j < m_M; ++j) {
                tmpd[2*j+0] = real(m[os+j]);
                tmpd[2*j+1] = imag(m[os+j]);
            }
            
			if (m_have_b0)
				NFFTTraits<NFFTType>::ITrafo
                    ((B0Plan&) m_b0_plan, (Solver&) m_solver, m_maxit, m_epsilon);
			else
				NFFTTraits<NFFTType>::ITrafo
                    ((Plan&)    m_plan,   (Solver&) m_solver, m_maxit, m_epsilon);
            
			tmpd = (NFFTRType*) m_solver.f_hat_iter;
			tmpt = (RT*) out.Ptr() + i*m_imgsz;

            std::copy (tmpd, tmpd+m_imgsz, tmpt);

			//TODO: b0 not 2D+1D+1D
			if (m_have_b0)
				for (size_t j = 0; j < out.Size(); ++j)
					out[j + i*m_imgsz / 2] *= std::polar<RT>((RT)1.,
							(RT)(-2. * PI * m_ts * m_b0[j] * m_w));

        }

/*        if (m_3rd_dim_cart && m_ncart > 1) { // Cartesian FT 3rd dim
        	int n = static_cast<int>(m_ncart);
        	out = permute (out, 2, 0, 1);
        	size_t cent = std::ceil(.5*m_ncart);
        	Matrix<T> tmp;
        	for (size_t i = 0; i < m_imgsz/2; ++i) { // fftshift colums
				tmp = Column (out,i);
				std::rotate(tmp.Begin(), tmp.Begin()+cent, tmp.End());
				Column (out, i, tmp);
        	}
        	CartPlan cp = FTTraits<T>::DFTPlanMany (1, &n, m_imgsz/2,
        			(FTType*)out.Ptr(),	(FTType*)out.Ptr(), FFTW_BACKWARD);
        	FTTraits<T>::Execute(cp);
        	FTTraits<T>::Destroy(cp);
        	for (size_t i = 0; i < m_imgsz/2; ++i) { // fftshift colums
				tmp = Column (out,i);
				std::rotate(tmp.Begin(), tmp.Begin()+cent, tmp.End());
				Column (out, i, tmp);
        	}
        	out = permute (out, 1, 2, 0)/sqrt((RT)m_ncart);
            }*/

        return out;
        
    }

    
    Plan& NFFTPlan() {return m_plan;}
    const Plan& NFFTPlan() const {return m_plan;}


    inline size_t Rank() const NOEXCEPT { return m_rank; }
    
    inline size_t ImageSize () const {return (size_t)m_plan.N[0];}
    inline size_t KSpaceSize () const {return (size_t)m_plan.M_total;}
    inline size_t Maxit () const {return m_maxit;}
    inline RT Alpha() const {return m_alpha;}
    inline RT Sigma() const {return m_sigma;}
    inline RT Epsilon() const {return m_epsilon;}

    virtual std::ostream& Print (std::ostream& os) const {
		Operator<T>::Print(os);
    	os << "    image size: rank(" << Rank() << ") side(" <<
            ImageSize() << ") nodes(" << KSpaceSize() << ")" << std::endl;
    	os << "    nfft: maxit(" << Maxit() << ") eps(" << Epsilon() <<	") alpha("
           << Alpha() << ") m("<< m_m <<") sigma(" << Sigma() << ")" << std::endl;
    	os << "    have_kspace(" << m_have_kspace << ") have_weights(" <<
            m_have_weights << ") have_b0(" << m_have_b0 << ")" << std::endl;
    	os << "    ft-threads(" << m_np << ")";
    	if (m_3rd_dim_cart)
    		os << " 3rd dimension (" << m_ncart << ") is Cartesian.";
    	return os;
    }

private:
    
    bool       m_initialised;   /**< @brief Memory allocated / Plans, well, planned! :)*/
    bool       m_have_pc, m_have_b0;
    
    size_t     m_rank;
    
    Matrix<T> m_pc;            /**< @brief Phase correction (applied after inverse trafo)*/
    Matrix<T> m_cpc;           /**< @brief Phase correction (applied before forward trafo)*/
    
    Matrix<RT>  m_b0;
    Matrix<RT>  m_t;

    RT m_min_t, m_max_t, m_min_b0, m_max_b0, m_w, m_ts;

    Vector<size_t> m_N;      /**< @brief Image matrix side length (incl. k_{\\omega})*/
    Vector<size_t> m_n;      /**< @brief Oversampling */
    Vector<RT> m_k, m_kw;

    size_t     m_M;             /**< @brief Number of k-space knots */
    size_t     m_maxit;         /**< @brief Number of Recon iterations (NFFT 3) */
    RT         m_epsilon;       /**< @brief Convergence criterium */
    RT         m_alpha, m_sigma;
    size_t     m_imgsz;
    
    Plan       m_plan;         /**< nfft  plan */
    B0Plan     m_b0_plan;
    Solver     m_solver;         /**< infft plan */
    CartPlan   m_cart_plan;
    
    bool       m_3rd_dim_cart, m_have_weights, m_have_kspace, m_per_slice_kspace;

    size_t     m_m, m_ncart;

    int m_np;

    //Window<RT>  m_win;

};


#endif




