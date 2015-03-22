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

#ifndef __NFFT_HPP__
#define __NFFT_HPP__

#include "NFFTTraits.hpp"
#include "FFTWTraits.hpp"
#include "Algos.hpp"
#include "Access.hpp"
#include "FT.hpp"
#include "CX.hpp"
#include "Creators.hpp"

/**
 * @brief Matrix templated ND non-equidistand Fourier transform with NFFT 3 (TU Chemnitz)<br/>
 *        T and single precision

 */
template <class T>
class NFFT : public FT<T> {

    typedef typename NFFTTraits<double>::Plan   Plan;
    typedef typename NFFTTraits<double>::B0Plan B0Plan;
    typedef typename NFFTTraits<double>::Solver Solver;
    typedef typename FTTraits<T>::Plan CartPlan;
	typedef typename FTTraits<T>::T FTType;

    typedef typename std::complex<T> CT;


public:

    /**
     * @brief         Default constructor
     */
    NFFT() NOEXCEPT :
        m_initialised (false),
        m_have_pc (false),
        m_imgsz (0),
        m_M (0),
        m_maxit (0),
        m_rank (0),
        m_m(0),
        m_have_b0(false),
        m_3rd_dim_cart(false),
        m_ncart(1) {};

    /**
     * @brief          Construct NFFT plans for forward and backward FT with credentials
     * 
     * @param  imsize  Matrix of side length of the image space
     * @param  nk      # k-space points
     * @param  m       Spatial cut-off of FT
     * @param  alpha   Oversampling factor
     * @param  b0      Off-resonance maps if available
     * @param  pc      Phase correction applied before forward or after adjoint transforms (default: empty)
     * @param  eps     Convergence criterium for inverse transform (default: 1.0e-7)
     * @param  maxit   Maximum # NFFT iterations (default: 3)
     */
    inline NFFT        (const Vector<size_t>& imsize, const size_t& nk, const size_t m = 1,
                        const T alpha = 1.0, const Matrix<T> b0 = Matrix<T>(1),
                        const Matrix< std::complex<T> > pc = Matrix< std::complex<T> >(1),
                        const T eps = 7.0e-4, const size_t maxit = 2) NOEXCEPT :
        m_have_b0(false), m_3rd_dim_cart(false), m_ncart(1){
        
        m_M     = nk;
        m_imgsz = 2;
        m_m = m;
        
        m_N = imsize;
        m_n = m_N;//ceil (alpha*m_N);
        
        m_rank = numel(imsize);
        
        m_imgsz = 2*prod(m_N);
        
        m_epsilon = eps;
        m_maxit   = maxit;
        
        NFFTTraits<double>::Init (m_N, m_M, m_n, m_m, m_plan, m_solver);
        
        if (pc.Size() > 1)
            m_have_pc = true;
        
        m_pc   = pc;
        m_cpc  = conj(pc);
        
        m_initialised = true;
        
    }
    
    
    inline NFFT (const Params& p) NOEXCEPT : m_have_b0(false), m_3rd_dim_cart(false),
        m_t (Matrix<T>(1)), m_b0 (Matrix<T>(1)), m_maxit(3), m_m(1), m_alpha(1.), m_epsilon(7.e-4f),
        m_sigma(1.0), m_ncart(1) {
                
        if (p.exists("nk")) {// Number of kspace samples
            try {
                m_M = unsigned_cast(p["nk"]);
            } catch (const boost::bad_any_cast& e) {
                printf ("**ERROR - NFFT: Numer of ksppace samples need to be specified\n%s\n", e.what());
                assert(false);
            }
        } else {
            printf ("**ERROR - NFFT: Numer of ksppace samples need to be specified\n");
            assert(false);
        }

        if (p.exists("3rd_dim_cart")) {
        	try {
        		m_3rd_dim_cart = p.Get<bool>("3rd_dim_cart");
        	} catch (const boost::bad_any_cast&) {
        		printf ("  WARNING - NFFT: Could not interpret input for Cartesian nature of 3rd dimension. \n");
        	}
        }

        if (p.exists("imsz")) {// Image domain size
            try {
                m_N = boost::any_cast<Vector<size_t> >(p["imsz"]);
            } catch (const boost::bad_any_cast& e) {
                printf ("**ERROR - NFFT: Image domain dimensions need to be specified\n%s\n", e.what());
                assert(false);
            }
        } else {
            printf ("**ERROR - NFFT: Image domain dimensions need to be specified\n");
            assert(false);
        }

        if (m_3rd_dim_cart) { // 3rd dimension is Cartesian
        	m_ncart = m_N.back();
        	m_N.PopBack();
        }
        m_n = m_N;
        
        if (p.exists("m")) {
            try {
                m_m = unsigned_cast (p["m"]);
            } catch (const boost::bad_any_cast&) {
                printf ("  WARNING - NFFT: Could not interpret input for oversampling factor m. "
                        "Defaulting to 1.\n");
            }
        }
        
        if (p.exists("alpha")) {
            try {
                m_alpha = fp_cast(p["alpha"]);
            } catch (const boost::bad_any_cast&) {
                printf ("  WARNING - NFFT: Could not interpret input for oversampling factor alpha. "
                        "Defaulting to 1.0\n");
            }
        }
        
        for (size_t i = 0; i < m_N.size(); ++i)
            m_n[i] = ceil(m_alpha*m_N[i]);
        
        m_rank  = m_N.size();
        m_imgsz = 2*prod(m_N);

        if (p.exists("epsilon")) {
            try {
                m_epsilon = fp_cast (p["epsilon"]);
            } catch (const boost::bad_any_cast&) {
                printf ("  WARNING - NFFT: Could not interpret input for convergence criterium epsilon. "
                        "Defaulting to 0.0007\n");
            }
        }
        
        if (p.exists("maxit")) {
            try {
                m_maxit = unsigned_cast (p["maxit"]);
            } catch (const boost::bad_any_cast&) {
                printf ("  WARNING - NFFT: Could not interpret input for maximum NFFT steps. "
                        "Defaulting to 3\n");
            }
        }

        if (p.exists("b0")) {
            try {
                m_b0 = p.Get<Matrix<T> >("b0");
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
                m_t = p.Get<Matrix<T> >("timing");
                if (numel(m_t) != m_M)
                    printf ("  WARNING - NFFT: kspace trajectory timing fit to trajectory. "
                            "Defaulting to flat b0 = 0\n");
                else
                    m_have_b0 = (m_have_b0 && true);
            } catch (const boost::bad_any_cast&) {
                printf ("  WARNING - NFFT: Could not interpret input for kspace trajectory timing. "
                        "Defaulting to flat b0 = 0\n");
            }
        }
        
        
        if (m_have_b0) {
            
            m_min_t  = min(m_t);
            m_max_t  = max(m_t);
            m_min_b0 = min(m_b0);
            m_max_b0 = max(m_b0);
            m_sigma  = 1.2;
            
            m_N.push_back(std::ceil(std::max(fabs(m_min_b0),fabs(m_max_b0)) *
                                    m_max_t-m_min_t/2.+(m_m)/(2.*m_sigma))*4.*m_sigma);
            
            if (m_N.back()%2!=0)
                m_N.back()++; // need even dimension

            m_n.push_back(m_N.back());
            m_n *= m_sigma;

            m_w = std::max(std::abs(m_min_b0),std::abs(m_max_b0))/(.5-((T) m_m)/m_N[2]);
            m_ts =  (m_min_t+m_max_t)/2.;
            T t    = ((m_max_t-m_min_t)/2.)/(.5-((T) (m_m))/m_N[2]);

            m_win = Window<T> (m_m, m_N.back(), m_sigma);
            
        	NFFTTraits<double>::Init (m_N, m_M, m_n, m_m, m_sigma, m_b0_plan, m_solver);

            for (size_t j = 0; j < m_N[0]*m_N[1]; ++j)
                m_b0_plan.w[j] = m_b0[j] / m_w;


        } else {

        	NFFTTraits<double>::Init (m_N, m_M, m_n, m_m, m_plan, m_solver);
        }
        
        if (p.exists("pc")) {
            try {
                m_pc = boost::any_cast<std::complex<T> >(p["pc"]);
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
    ~NFFT () NOEXCEPT {
        if (m_initialised)
        	if (m_have_b0)
        		NFFTTraits<double>::Finalize (m_b0_plan, m_solver);
        	else
        		NFFTTraits<double>::Finalize (m_plan, m_solver);
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
        m_win         = ft.m_win;
        m_sigma       = ft.m_sigma;
        m_3rd_dim_cart = ft.m_3rd_dim_cart;
        m_ncart       = ft.m_ncart;
        if (m_have_b0)
        	NFFTTraits<double>::Init (m_N, m_M, m_n, m_m, m_sigma, m_b0_plan, m_solver);
        else
        	NFFTTraits<double>::Init (m_N, m_M, m_n, m_m, m_plan, m_solver);
        return *this;
        
    }
    
    /**
     * @brief      Assign k-space 
     * 
     * @param  k   Kspace trajectory
     */
    inline void 
    KSpace (const Matrix<T>& k) NOEXCEPT {
        if (m_have_b0) { // +1D for omega
            for (size_t j = 0; j < m_b0_plan.M_total; ++j) {
                m_b0_plan.plan.x[3*j+0] = real(k[2*j+0]);
                m_b0_plan.plan.x[3*j+1] = imag(k[2*j+1]);
                m_b0_plan.plan.x[3*j+2] = (m_t[j]-m_ts)*m_w/m_N.back();
            }
        } else {
            assert (k.Size() == m_plan.M_total*m_rank);
            std::copy (k.Begin(), k.End(), m_plan.x);
        }
    }
    
    
    /**
     * @brief      Assign k-space weigths (jacobian of trajectory with regards to time) 
     * 
     * @param  w   Weights
     */
    inline void 
    Weights (const Matrix<T>& w) NOEXCEPT {

    	if (m_have_b0)
    		assert (w.Size() == m_b0_plan.M_total);
    	else
    		assert (w.Size() == m_plan.M_total);

        std::copy (w.Begin(), w.End(), m_solver.w);

        if (m_have_b0) {
            NFFTTraits<double>::Weights (m_b0_plan.plan, m_solver, m_rank);
            NFFTTraits<double>::Psi (m_b0_plan.plan);
        } else {
        	NFFTTraits<double>::Weights (m_plan, m_solver, m_rank);
        	NFFTTraits<double>::Psi (m_plan);
        }

    }
    
    
    /**
     * @brief    Forward transform
     *
     * @param  m To transform
     * @return   Transform
     */
    Matrix< std::complex<T> >
    Trafo       (const Matrix< std::complex<T> >& m) const NOEXCEPT {

		double* tmpd;
		T* tmpt;
        Matrix< std::complex<T> > out (m_M * ((m_3rd_dim_cart && m_ncart > 1) ? m_ncart : 1), 1);
        Matrix< std::complex<T> > tmpm = m;

        if (m_3rd_dim_cart && m_ncart > 1) { // Cartesian FT 3rd dim
        	int n = static_cast<int>(m_ncart);
        	tmpm = permute (tmpm, 2, 0, 1);
        	size_t cent = floor(m_ncart/2);
        	Matrix< std::complex<T> > tmp;
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
        	tmpm = permute (tmpm, 1, 2, 0);
        }


        for (size_t i = 0; i < m_ncart; ++i) {

			tmpd = (double*) m_plan.f_hat;
			tmpt = (T*) tmpm.Ptr() + i*m_imgsz;

			//TODO: b0 not 2D+1D+1D
			if (m_have_b0)
				for (size_t j = 0; j < m.Size(); ++j) {
					CT val = tmpm[j] * std::polar<T> ((T)1., (T)(2. * PI * m_ts * m_b0[j] * m_w));
					tmpd[2*j+0] = (T)real(val);
					tmpd[2*j+1] = (T)imag(val);
				}
			else
				std::copy (tmpt, tmpt+m_imgsz, tmpd);

			if (m_have_b0)
				NFFTTraits<double>::Trafo (m_b0_plan);
			else
				NFFTTraits<double>::Trafo (m_plan);

			tmpd = (double*) m_plan.f;
			tmpt = (T*) out.Ptr() + i*2*m_M;
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
    Matrix< std::complex<T> >
    Adjoint     (const Matrix< std::complex<T> >& m) const NOEXCEPT {

        Vector<size_t> N = m_N;
        double* tmpd;
        T* tmpt;

        if (m_have_b0)
            N.PopBack();
        if (m_3rd_dim_cart && m_ncart > 1) // Cartesian FT 3rd dim
        	N.PushBack(m_ncart);

        Matrix<std::complex<T> > out (N);
        for (size_t i = 0; i < m_ncart; ++i) {

			tmpd = (double*) m_solver.y;
			tmpt = (T*) m.Ptr() + i*2*m_M;
			std::copy (tmpt, tmpt+2*m_M, tmpd);

			if (m_have_b0)
				NFFTTraits<double>::ITrafo ((B0Plan&) m_b0_plan, (Solver&) m_solver, m_maxit, m_epsilon);
			else
				NFFTTraits<double>::ITrafo (  (Plan&)    m_plan, (Solver&) m_solver, m_maxit, m_epsilon);

			tmpd = (double*) m_solver.f_hat_iter;
			tmpt = (T*) out.Ptr() + i*m_imgsz;
			std::copy (tmpd, tmpd+m_imgsz, tmpt);

			//TODO: b0 not 2D+1D+1D
			if (m_have_b0)
				for (size_t j = 0; j < out.Size(); ++j)
					out[j + i*m_imgsz / 2] *= std::polar<T>((T)1., (T)(-2. * PI * m_ts * m_b0[j] * m_w));

        }

        if (m_3rd_dim_cart && m_ncart > 1) { // Cartesian FT 3rd dim
        	int n = static_cast<int>(m_ncart);
        	out = permute (out, 2, 0, 1);
        	size_t cent = floor(m_ncart/2);
        	Matrix< std::complex<T> > tmp;
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
        	out = permute (out, 1, 2, 0);
        }

        return out;
        
    }

    

    inline size_t Rank() const NOEXCEPT { return m_rank; }
    
private:
    
    bool       m_initialised;   /**< @brief Memory allocated / Plans, well, planned! :)*/
    bool       m_have_pc, m_have_b0;
    
    size_t     m_rank;
    
    Matrix<CT> m_pc;            /**< @brief Phase correction (applied after inverse trafo)*/
    Matrix<CT> m_cpc;           /**< @brief Phase correction (applied before forward trafo)*/
    
    Matrix<T>  m_b0;
    Matrix<T>  m_t;

    T m_min_t, m_max_t, m_min_b0, m_max_b0, m_w, m_ts;

    Vector<size_t> m_N;      /**< @brief Image matrix side length (incl. k_{\\omega})*/
    Vector<size_t> m_n;      /**< @brief Oversampling */

    size_t     m_M;             /**< @brief Number of k-space knots */
    size_t     m_maxit;         /**< @brief Number of Recon iterations (NFFT 3) */
    T          m_epsilon;       /**< @brief Convergence criterium */
    T          m_alpha, m_sigma;
    size_t     m_imgsz;
    
    Plan       m_plan;         /**< nfft  plan */
    B0Plan     m_b0_plan;
    Solver     m_solver;         /**< infft plan */
    CartPlan   m_cart_plan;
    
    bool       m_3rd_dim_cart;

    size_t     m_m, m_ncart;

    Window<T>  m_win;

};


#endif




