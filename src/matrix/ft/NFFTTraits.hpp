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

#include <math.h>
#include <stdlib.h>
#include "Complex.hpp"

#include "nfft3util.h"
#include "nfft3.h"

#include "Container.hpp"

#define USE_NFFT_32_NAMING 1

#ifndef USE_NFFT_32_NAMING
    #define nfft_mv_plan_complex mv_plan_complex
    #define nfftf_mv_plan_complex mv_plan_complex
#endif

template <class T>
struct NFFTTraits { };


#ifdef USE_NFFT_32_NAMING

template <>
struct NFFTTraits<float> {

	typedef nfftf_plan           Plan;    /**< @brief nfft plan (float precision) */
	typedef solverf_plan_complex Solver;  /**< @brief nfft solver plan (float precision) */

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
	Init  (const int d, int* N, const int M, int* n, const int m, nfftf_plan& np, solverf_plan_complex& inp) {

		fftwf_init_threads();
		
#ifdef _OPENMP
		fftwf_import_wisdom_from_filename("codeare_threads.plan");
#else
		fftwf_import_wisdom_from_filename("codeare_single.plan");
#endif

		nfftf_init_guru
			(&np, d, N, M, n, m,
					NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT |
					PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
					FFTW_MEASURE| FFTW_DESTROY_INPUT);


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
	ITrafo              (nfftf_plan& np, solverf_plan_complex& spc, const int maxiter = 3, const double epsilon = 3e-7) {
		
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
	Trafo                (const nfftf_plan& np) {
		
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
	Adjoint              (const nfftf_plan& np) {
		
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
	Weights              (const nfftf_plan& np, const solverf_plan_complex& spc) {
		
		int j, k, z, N = np.N[0];
		
		if (spc.flags & PRECOMPUTE_DAMP) {
			if (np.d == 3) {
				for (j = 0; j < N; j++) {
					float    j2 = j - N/2; 
					for (k = 0; k < N; k++) {
						float    k2 = k - N/2;
						for (z = 0; z < N; z++) {
							float    z2 = z - N/2;
							float r  = sqrt(j2*j2+k2*k2+z2*z2);
							spc.w_hat[z*N*N+j*N+k] = (r > (float) N/2) ? 0.0 : 1.0;
						}
					}
				}
			} else {
				for (j = 0; j < N; j++) {
					float    j2 = j-N/2;
					for (k = 0; k < N; k++) {
						float    k2 = k-N/2;
						double r  = sqrt(j2*j2+k2*k2);
						spc.w_hat[j*N+k]       = (r > (float) N/2) ? 0.0 : 1.0;
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
	Psi                  (nfftf_plan& np) {
		
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
	Finalize             (nfftf_plan& np, solverf_plan_complex& spc) {
		
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
struct NFFTTraits<double> {

	typedef nfft_plan           Plan;    /**< @brief nfft plan (float precision) */
	typedef solver_plan_complex Solver;  /**< @brief nfft solver plan (float precision) */

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
	Init  (const container<int>& N, size_t M, const container<int>& n,
			int m, nfft_plan& np, solver_plan_complex& inp) {


		container<int> _N(N), _n(n);
		int _d (N.size()), _M(M), _m(m);


#ifdef _OPENMP
		fftw_import_wisdom_from_filename("codeare_threads.plan");
#else
		fftw_import_wisdom_from_filename("codeare_single.plan");
#endif

		nfft_init_guru 
			(&np, _d, _N.ptr(), _M, _n.ptr(), _m,
			 NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT | 
             PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
			 FFTW_MEASURE| FFTW_DESTROY_INPUT);
		
		solver_init_advanced_complex 
			(&inp, (nfft_mv_plan_complex*) &np, CGNR | PRECOMPUTE_DAMP | PRECOMPUTE_WEIGHT);
		
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
	ITrafo              (nfft_plan& np, solver_plan_complex& spc, size_t maxiter = 3, double epsilon = 3e-7) {
		
		int l;
		
		/* init some guess */
        std::fill_n ((double*)spc.f_hat_iter, 2*np.N_total, 0.);
		
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
	Trafo                (const nfft_plan& np) {
		
		nfft_trafo ((nfft_plan*) &np);
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
	Adjoint              (const nfft_plan& np) {
		
		nfft_adjoint ((nfft_plan*) &np);
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
	Weights              (const nfft_plan& np, const solver_plan_complex& spc) {
		
		int j, k, z, N = np.N[0], N2 = N*N, NH = .5*N;
		float k2, j2, z2;
		
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
	Psi                  (nfft_plan& np) {
		
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
	Finalize             (nfft_plan& np, solver_plan_complex& spc) {
		
		solver_finalize_complex(&spc);
		nfft_finalize(&np);
		return 0;
		
	}
	
	
};



#endif //__NFFT_TRAITS_HPP__
