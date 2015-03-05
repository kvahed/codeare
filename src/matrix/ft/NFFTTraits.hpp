/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
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

#ifndef __NFFT_TRAITS_HPP__
#define __NFFT_TRAITS_HPP__

#include <cmath>
#include <math.h>
#include <stdlib.h>
#include "Complex.hpp"

#include "nfft3util.h"
#include "nfft3.h"

#include "Vector.hpp"

#define USE_NFFT_32_NAMING 1

#ifndef USE_NFFT_32_NAMING
    #define nfft_mv_plan_complex mv_plan_complex
    #define nfftf_mv_plan_complex mv_plan_complex
#endif

enum nfft_startegy {WITHOUT_B0, WITH_B0};

template <class T, nfft_startegy S = WITHOUT_B0>
struct NFFTTraits { };

#ifdef USE_NFFT_32_NAMING


template <>
struct NFFTTraits<float, WITHOUT_B0> {

	typedef nfftf_plan           Plan;    /**< @brief nfft plan (float precision) */
	typedef solverf_plan_complex Solver;  /**< @brief nfft solver plan (float precision) */
	typedef float                T;

	/**
	 * @brief            Initialise plan
	 *
	 * @param  d         Number of dimension (i.e. {1..3})
	 * @param  N         Actual dimensions 
	 * @param  M         Number of k-space samples
	 * @param  n         Oversampled N
	 * @param  m         Spatial cutoff
	 * @param  np        Forward FT plan
	 * @param  inp       Inverse FT plan
	 *
	 * @return success
	 */
	inline static int
	Init  (const int d, int* N, const int M, int* n, const int m, nfftf_plan& np,
			solverf_plan_complex& inp) NOEXCEPT {

		fftwf_init_threads();
		
#ifdef _OPENMP
		fftwf_import_wisdom_from_filename("codeare_threads.plan");
#else
		fftwf_import_wisdom_from_filename("codeare_single.plan");
#endif

		nfftf_init_guru
			(&np, d, N, M, n, m,
					NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT | PRE_PHI_HUT |
					PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT |
					FFT_OUT_OF_PLACE, FFTW_MEASURE| FFTW_DESTROY_INPUT);


		solverf_init_advanced_complex
			(&inp, (nfftf_mv_plan_complex*) &np, CGNR | PRECOMPUTE_DAMP | PRECOMPUTE_WEIGHT);

#ifdef _OPENMP
		fftwf_export_wisdom_to_filename("codeare_threads.plan");
#else
		fftwf_export_wisdom_to_filename("codeare_single.plan");
#endif

		return 0;
		
	}
    
    

	/**
	 * @brief            Inverse FT
	 * 
	 * @param  np        NFFT plan
	 * @param  spc       iNFFT plan
	 * @param  maxiter   Maximum NFFT iterations        
	 * @param  epsilon   Convergence criterium
	 *
	 * @return           Success
	 */
	inline static int
	ITrafo              (nfftf_plan& np, solverf_plan_complex& spc, const int maxiter = 3,
			const T epsilon = 3e-7) NOEXCEPT {
		
		int k, l;
		
		/* init some guess */
		for (k = 0; k < np.N_total; k++) {
			spc.f_hat_iter[k][0] = 0.0;
			spc.f_hat_iter[k][1] = 0.0;
		}
		
		/* inverse trafo */
		solverf_before_loop_complex(&spc);
		
		for (l = 0; l < maxiter; l++) {
			if (spc.dot_r_iter < epsilon) 
				break;
			solverf_loop_one_step_complex(&spc);
		}
		
		return 0;
		
	}
	


	/**
	 * @brief            Forward FT
	 *
	 * @param  np        NFFT plan
	 *
	 * @return           Success
	 */
	inline static int
	Trafo                (const nfftf_plan& np) NOEXCEPT {
		
		nfftf_trafo ((nfftf_plan*) &np);
		return 0;
		
	}
	
	

	/**
	 * @brief            Adjoint FT
	 *
	 * @param  np        NFFT plan
	 *
	 * @return           Success
	 */
	inline static int
	Adjoint              (const nfftf_plan& np) NOEXCEPT {
		
		nfftf_adjoint ((nfftf_plan*) &np);
		return 0;
		
	}


	/**
	 * @brief            Set weights
	 *
	 * @param  np        Plan
	 * @param  spc       Solver plan
	 * @return           Success
	 */
	inline static int
	Weights              (const nfftf_plan& np, const solverf_plan_complex& spc) NOEXCEPT {

		int j, k, z, N = np.N[0];

		if (spc.flags & PRECOMPUTE_DAMP) {
			if (np.d == 3) {
				for (j = 0; j < N; j++) {
					T j2 = j - N/2;
					for (k = 0; k < N; k++) {
						T k2 = k - N/2;
						for (z = 0; z < N; z++) {
							T z2 = z - N/2;
							T r  = sqrt(j2*j2+k2*k2+z2*z2);
							spc.w_hat[z*N*N+j*N+k] = (r > (T) N/2) ? 0.0 : 1.0;
						}
					}
				}
			} else {
				for (j = 0; j < N; j++) {
					T    j2 = j-N/2;
					for (k = 0; k < N; k++) {
						T    k2 = k-N/2;
						T r  = sqrt(j2*j2+k2*k2);
						spc.w_hat[j*N+k]       = (r > (T) N/2) ? 0.0 : 1.0;
					}
				}
			}
		}

		return np.M_total;

	}



	/**
	 * @brief            Precompute PSI
	 *
	 * @param  np        NFFT plan
	 * @return           Success
	 */
	inline static int
	Psi                  (nfftf_plan& np) NOEXCEPT {
		
		/* precompute full psi */
		if(np.nfft_flags & PRE_PSI)
			nfftf_precompute_one_psi(&np);
		
		/* precompute full psi */
		if(np.nfft_flags & PRE_FULL_PSI)
			nfftf_precompute_full_psi(&np);

		return 0;
		
	}


	
	/**
	 * @brief            Finalise plans
	 *
	 * @param  np        Plan
	 * @param  spc       Solver plan
	 * @return           Success
	 */
	inline static int
	Finalize             (nfftf_plan& np, solverf_plan_complex& spc) NOEXCEPT {
		
		solverf_finalize_complex(&spc);
		nfftf_finalize(&np);
#ifdef _OPENMP
		fftwf_cleanup_threads();
#endif
		return 0;
		
	}
	
	
};

#endif


template <>
struct NFFTTraits<double, WITHOUT_B0> {

	typedef nfft_plan           Plan;    /**< @brief nfft plan (double precision) */
	typedef solver_plan_complex Solver;  /**< @brief nfft solver plan (double precision) */
	typedef double              T;

	/**
	 * @brief            Initialise plan
	 *
	 * @param  d         Number of dimension (i.e. {1..3})
	 * @param  N         Actual dimensions
	 * @param  M         Number of k-space samples
	 * @param  n         Oversampled N
	 * @param  m         Spatial cutoff
	 * @param  np        Forward FT plan
	 * @param  inp       Inverse FT plan
	 *
	 * @return success
	 */
	inline static int
	Init  (const Vector<int>& N, size_t M, const Vector<int>& n, int m, Plan& np, Solver& inp) NOEXCEPT {


		Vector<int> _N(N), _n(n);
		int _d (N.size()), _M(M), _m(m);
		unsigned nfft_flags, infft_flags, solver_flags;

		infft_flags  = FFTW_MEASURE| FFTW_DESTROY_INPUT;
		solver_flags = CGNR | PRECOMPUTE_DAMP | PRECOMPUTE_WEIGHT;
		nfft_flags   = NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT | PRE_PHI_HUT |
				PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE;

#ifdef _OPENMP
		fftw_import_wisdom_from_filename("codeare_threads.plan");
#else
		fftw_import_wisdom_from_filename("codeare_single.plan");
#endif

		nfft_init_guru (&np, _d, _N.ptr(), _M, _n.ptr(), _m, nfft_flags, infft_flags);
		solver_init_advanced_complex (&inp, (nfft_mv_plan_complex*) &np, solver_flags);

#ifdef _OPENMP
		fftw_export_wisdom_to_filename("codeare_threads.plan");
#else
		fftw_export_wisdom_to_filename("codeare_single.plan");
#endif

		return 0;

	}

	/**
	 * @brief            Initialise plan
	 *
	 * @param  d         Number of dimension (i.e. {1..3})
	 * @param  N         Actual dimensions
	 * @param  M         Number of k-space samples
	 * @param  n         Oversampled N
	 * @param  m         Spatial cutoff
	 * @param  np        Forward FT plan
	 * @param  inp       Inverse FT plan
	 *
	 * @return success
	 */
	inline static int
	Init  (const Vector<int>& N, size_t M, const Vector<int>& n, int m, const Vector<T>& b0,
			const Vector<T>& t, Plan& np, Solver& inp) NOEXCEPT {


		Vector<int> _N(N), _n(n);
		int _d (N.size()), _M(M), _m(m);
		unsigned nfft_flags, infft_flags, solver_flags;

		infft_flags  = FFTW_MEASURE| FFTW_DESTROY_INPUT;
		solver_flags = CGNR | PRECOMPUTE_DAMP | PRECOMPUTE_WEIGHT;
		nfft_flags   = NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT | PRE_PHI_HUT |
				PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE;

/*
		  ths->N3=N[2];
		  ths->sigma3=sigma;
		  nfft_init_guru(&ths->plan,3,N,M,n,m,nfft_flags,fftw_flags);
		  ths->N_total = N[0]*N[1];
		  ths->M_total = ths->plan.M_total;
		  ths->f = ths->plan.f;
		  ths->f_hat = (double _Complex*) nfft_malloc(ths->N_total*sizeof(double _Complex));
		  ths->w = (double*) nfft_malloc(ths->N_total*sizeof(double));

		  ths->mv_trafo = (void (*) (void* ))mri_inh_3d_trafo;
		  ths->mv_adjoint = (void (*) (void* ))mri_inh_3d_adjoint;
*/

#ifdef _OPENMP
		fftw_import_wisdom_from_filename("codeare_threads.plan");
#else
		fftw_import_wisdom_from_filename("codeare_single.plan");
#endif

		nfft_init_guru (&np, _d, _N.ptr(), _M, _n.ptr(), _m, nfft_flags, infft_flags);
		solver_init_advanced_complex (&inp, (nfft_mv_plan_complex*) &np, solver_flags);

#ifdef _OPENMP
		fftw_export_wisdom_to_filename("codeare_threads.plan");
#else
		fftw_export_wisdom_to_filename("codeare_single.plan");
#endif

		return 0;

	}



	/**
	 * @brief            Inverse FT
	 *
	 * @param  np        NFFT plan
	 * @param  spc       iNFFT plan
	 * @param  maxiter   Maximum NFFT iterations
	 * @param  epsilon   Convergence criterium
	 *
	 * @return           Success
	 */
	inline static int
	ITrafo              (Plan& np, Solver& spc, size_t maxiter = 3,	T epsilon = 3e-7) NOEXCEPT {

		size_t l;

		/* init some guess */
        std::fill_n ((T*)spc.f_hat_iter, 2*np.N_total, 0.);

		/* inverse trafo */
		solver_before_loop_complex(&spc);

		for (l = 0; l < maxiter; l++) {
			if (spc.dot_r_iter < epsilon)
				break;
			solver_loop_one_step_complex(&spc);
		}

		return 0;

	}



	/**
	 * @brief            Forward FT
	 *
	 * @param  np        NFFT plan
	 *
	 * @return           Success
	 */
	inline static int
	Trafo                (const Plan& np) NOEXCEPT {

		nfft_trafo ((Plan*) &np);
		return 0;

	}



	/**
	 * @brief            Adjoint FT
	 *
	 * @param  np        NFFT plan
	 *
	 * @return           Success
	 */
	inline static int
	Adjoint              (const Plan& np) NOEXCEPT {

		nfft_adjoint ((Plan*) &np);
		return 0;

	}


	/**
	 * @brief            Set weights
	 *
	 * @param  np        Plan
	 * @param  spc       Solver plan
	 * @return           Success
	 */
	inline static int
	Weights              (const Plan& np, const Solver& spc) NOEXCEPT {

		int j, k, z, N = np.N[0], N2 = N*N, NH = .5*N;
		T k2, j2, z2;

		if (spc.flags & PRECOMPUTE_DAMP)
			if (np.d == 3) {
				for (j = 0; j < N; ++j) {
					j2 = j - NH;
					j2 *= j2;
					for (k = 0; k < N; ++k) {
						k2 = k - NH;
						k2 *= k2;
						for (z = 0; z < N; ++z) {
							z2 = z - NH;
							spc.w_hat[z*N2+j*N+k] = (sqrt(j2+k2+z2*z2) > NH) ? 0. : 1.;
						}
					}
				}
			} else {
				for (j = 0; j < N; j++) {
					j2 = j-NH;
					j2 *= j2;
					for (k = 0; k < N; k++) {
						k2 = k-NH;
						spc.w_hat[j*N+k] = (sqrt(j2+k2*k2) > NH) ? 0. : 1.;
					}
				}
			}

		return np.M_total;

	}



	/**
	 * @brief            Precompute PSI
	 *
	 * @param  np        NFFT plan
	 * @return           Success
	 */
	inline static int
	Psi                  (Plan& np) NOEXCEPT {

		/* precompute full psi */
		if(np.nfft_flags & PRE_PSI)
			nfft_precompute_one_psi(&np);

		/* precompute full psi */
		if(np.nfft_flags & PRE_FULL_PSI)
			nfft_precompute_full_psi(&np);

		return 0;

	}



	/**
	 * @brief            Finalise plans
	 *
	 * @param  np        Plan
	 * @param  spc       Solver plan
	 * @return           Success
	 */
	inline static int
	Finalize             (Plan& np, Solver& spc) NOEXCEPT {

		solver_finalize_complex(&spc);
		nfft_finalize(&np);
		return 0;

	}


};


#ifdef FALSE
template <>
struct NFFTTraits<double, WITH_B0> {

	typedef mri_inh_2d1d_plan     Plan;    /**< @brief nfft plan (double precision) */
	typedef solver_plan_complex Solver;  /**< @brief nfft solver plan (double precision) */
	typedef double              T;

	/**
	 * @brief            Initialise plan
	 *
	 * @param  d         Number of dimension (i.e. {1..3})
	 * @param  N         Actual dimensions
	 * @param  M         Number of k-space samples
	 * @param  n         Oversampled N
	 * @param  m         Spatial cutoff
	 * @param  np        Forward FT plan
	 * @param  inp       Inverse FT plan
	 *
	 * @return success
	 */
	inline static int
	Init  (const Vector<int>& N, size_t M, const Vector<int>& n, int m,	const Vector<T>& time,
			const Vector<T>& b0, Plan& np, Solver& inp) NOEXCEPT {


		Vector<int> _N(N), _n(n);
		int _d (N.size()), _M(M), _m(m);
		T sigma = 1.25, min_time, max_time, min_inh, max_inh, t, w, ts;

		unsigned nfft_flags, infft_flags, solver_flags;

		infft_flags  = FFTW_MEASURE| FFTW_DESTROY_INPUT;
		solver_flags = CGNR | PRECOMPUTE_DAMP | PRECOMPUTE_WEIGHT;
		nfft_flags   = NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT | PRE_PHI_HUT |
				PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE;

		min_time = INT_MAX;
		max_time = INT_MIN;
		for (size_t j=0; j < M; ++j) {
		    if (time[j] < min_time) min_time = time[j];
		    if (time[j] > max_time) max_time = time[j];
		}

		min_inh = INT_MAX;
		max_inh = INT_MIN;
		for (size_t j=0; j < _N[0]*_N[1]; ++j) {
		    if (b0[j] < min_inh) min_inh = b0[j];
		    if (b0[j] > max_inh) max_inh = b0[j];
		}

		_N.push_back(std::ceil(std::max(fabs(min_inh),fabs(max_inh)) *
				max_time-min_time/2.+(m)/(2.*sigma))*4.*sigma);

		if(_N[2]%2!=0)
		  _N[2]++;

		_n.push_back(_N[2]);

		ts = (min_time+max_time)/2.;
		t  = ((max_time-min_time)/2.0)/(0.5-((T) (m))/_N[2]);
		w  = _N[2]/t;

		for (size_t j=0; j < _N[0]*_N[1]; ++j)
			np.w[j]= b0[j]/w;

		for (size_t j = 0; j < np.M_total; ++j)
			np.t[j] = (time[j]-ts)/t;

#ifdef _OPENMP
		fftw_import_wisdom_from_filename("codeare_threads.plan");
#else
		fftw_import_wisdom_from_filename("codeare_single.plan");
#endif

		mri_inh_2d1d_init_guru (&np, _N.ptr(), _M, _n.ptr(), _m, sigma, nfft_flags, infft_flags);
		solver_init_advanced_complex (&inp, (nfft_mv_plan_complex*) &np, solver_flags);

#ifdef _OPENMP
		fftw_export_wisdom_to_filename("codeare_threads.plan");
#else
		fftw_export_wisdom_to_filename("codeare_single.plan");
#endif

		return 0;

	}


	/**
	 * @brief            Inverse FT
	 *
	 * @param  np        NFFT plan
	 * @param  spc       iNFFT plan
	 * @param  maxiter   Maximum NFFT iterations
	 * @param  epsilon   Convergence criterium
	 *
	 * @return           Success
	 */
	inline static int
	ITrafo              (Plan& np, Solver& spc, size_t maxiter = 3,	T epsilon = 3e-7) NOEXCEPT {

		size_t l;

		/* init some guess */
        std::fill_n ((T*)spc.f_hat_iter, 2*np.N_total, 0.);

		/* inverse trafo */
		solver_before_loop_complex(&spc);

		for (l = 0; l < maxiter; l++) {
			if (spc.dot_r_iter < epsilon)
				break;
			solver_loop_one_step_complex(&spc);
		}

		//for (size_t j = 0; j < 2*)
		//spc.f_hat_iter[j]*=cexp(-2.0*_Complex_I*PI*Ts*my_plan.w[j]*W);

		return 0;

	}


	/**
	 * @brief            Forward FT
	 *
	 * @param  np        NFFT plan
	 *
	 * @return           Success
	 */
	inline static int
	Trafo                (const Plan& np) NOEXCEPT {

		T ts =  .5 * (np.t[0]+np.t[np.M_total-1]);
		T t  = (.5 * (np.t[0]-np.t[np.M_total-1]))/(0.5-((T) np.plan.m)/np.plan.N[2]);
		T w  = np.plan.N[2]/t;

		for(size_t j = 0; j < np.plan.N[0]*np.plan.N[1]; ++j)
		    np.f_hat[j] *= std::polar<T> (1., 2.*PI*ts*np.w[j]/w);

		mri_inh_2d1d_trafo ((Plan*) &np);

		return 0;

	}



	/**
	 * @brief            Adjoint FT
	 *
	 * @param  np        NFFT plan
	 *
	 * @return           Success
	 */
	inline static int
	Adjoint              (const Plan& np) NOEXCEPT {

		mri_inh_2d1d_adjoint ((Plan*) &np);
		return 0;

	}


	/**
	 * @brief            Set weights
	 *
	 * @param  np        Plan
	 * @param  spc       Solver plan
	 * @return           Success
	 */
	inline static int
	Weights              (const Plan& np, const Solver& spc) NOEXCEPT {

		int j, k, z, N = np.plan.N[0], N2 = N*N, NH = .5*N;
		T k2, j2, z2;

		if (spc.flags & PRECOMPUTE_DAMP)
			if (np.plan.d == 3) {
				for (j = 0; j < N; ++j) {
					j2 = j - NH;
					j2 *= j2;
					for (k = 0; k < N; ++k) {
						k2 = k - NH;
						k2 *= k2;
						for (z = 0; z < N; ++z) {
							z2 = z - NH;
							spc.w_hat[z*N2+j*N+k] = (sqrt(j2+k2+z2*z2) > NH) ? 0. : 1.;
						}
					}
				}
			} else {
				for (j = 0; j < N; j++) {
					j2 = j-NH;
					j2 *= j2;
					for (k = 0; k < N; k++) {
						k2 = k-NH;
						spc.w_hat[j*N+k] = (sqrt(j2+k2*k2) > NH) ? 0. : 1.;
					}
				}
			}

		return np.M_total;

	}


	/**
	 * @brief            Precompute PSI
	 *
	 * @param  np        NFFT plan
	 * @return           Success
	 */
	inline static int
	Psi                  (Plan& np) NOEXCEPT {

		/* precompute full psi */
		if(np.plan.nfft_flags & PRE_PSI)
			nfft_precompute_one_psi(&np.plan);

		/* precompute full psi */
		if(np.plan.nfft_flags & PRE_FULL_PSI)
			nfft_precompute_full_psi(&np.plan);

		return 0;

	}

	/**
	 * @brief            Finalise plans
	 *
	 * @param  np        Plan
	 * @param  spc       Solver plan
	 * @return           Success
	 */
	inline static int
	Finalize             (Plan& np, Solver& spc) NOEXCEPT {

		solver_finalize_complex(&spc);
		mri_inh_2d1d_finalize(&np);
		return 0;

	}


};
#endif

#endif //__NFFT_TRAITS_HPP__
