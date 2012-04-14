/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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

#ifndef __NFFT_STUB_HPP__
#define __NFFT_STUB_HPP__

#include <math.h>
#include <stdlib.h>
#include "Complex.hpp"

#include "nfft3util.h"
#include "nfft3.h"

/**
 * @brief Interface to <a href="http://www-user.tu-chemnitz.de/~potts/nfft/" target="NFFT">NFFT 3 library</a>
 */
namespace nnfft {

		
	/**
	 * @brief            Initialise plan
	 *
	 * @param  d         Number of dimension (i.e. {1..3})
	 * @param  N         Actual dimensions 
	 * @param  M         Number of k-space samples
	 * @param  n         Oversampled N
	 * @param  m         Spatial cutoff
	 * @param  fnp       Forward FT plan
	 * @param  inp       Inverse FT plan
	 *
	 * @return success
	 */
	inline static int
	init                 (const int d, int* N, const int M, int* n, const int m, nfft_plan& np, solver_plan_complex& inp) {
		
		nfft_init_guru 
			(&np, d, N, M, n, m,
			 PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
			 FFTW_MEASURE| FFTW_DESTROY_INPUT);
		
		solver_init_advanced_complex 
			(&inp, (mv_plan_complex*)&np, CGNR | PRECOMPUTE_DAMP | PRECOMPUTE_WEIGHT);
		
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
	ift                  (nfft_plan& np, solver_plan_complex& spc, const int maxiter = 3, const double epsilon = 3e-7) {
		
		int k, l;
		
		/* init some guess */
		for (k = 0; k < np.N_total; k++) {
			spc.f_hat_iter[k][0] = 0.0;
			spc.f_hat_iter[k][1] = 0.0;
		}
		
		//double t = nfft_second();
		
		/* inverse trafo */
		solver_before_loop_complex(&spc);
		
		for (l = 0; l < maxiter; l++) {
			
			if (spc.dot_r_iter < epsilon) 
				break;
			
			solver_loop_one_step_complex(&spc);
			
#ifdef NVERBOSE		
			fprintf( stderr, "%e,  %i of %i\n", sqrt(spc.dot_r_iter), l+1, maxiter);
#endif		
			
		}
		
		//t = nfft_second()-t;
		
#ifdef NVERBOSE
		//fprintf(stderr,"nfft time: %.4f seconds.\n",t);
#endif
		
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
	ft                   (const nfft_plan& np) {
		
		nfft_trafo ((nfft_plan*) &np);
		return 0;
		
	}
	
	
	/**
	 * @brief            Adjoint FT
	 *
	 * @param  np        NFFT plan
	 
	 * @return           Success
	 */
	inline static int
	adjoint              (const nfft_plan& np) {
		
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
	int
	weights              (const nfft_plan& np, const solver_plan_complex& spc) {
		
		int j, k, z, N = np.N[0];
		
		if (spc.flags & PRECOMPUTE_DAMP) {
			if (np.d == 3) {
				for (j = 0; j < N; j++) {
					int    j2 = j - N/2; 
					for (k = 0; k < N; k++) {
						int    k2 = k - N/2;
						for (z = 0; z < N; z++) {
							int    z2 = z - N/2;
							double r  = sqrt(j2*j2+k2*k2+z2*z2);
							spc.w_hat[z*N*N+j*N+k] = (r > (double) N/2) ? 0.0 : 1.0;
						}
					}
				}
			} else {
				for (j = 0; j < N; j++) {
					int    j2 = j-N/2;
					for (k = 0; k < N; k++) {
						int    k2 = k-N/2;
						double r  = sqrt(j2*j2+k2*k2);
						spc.w_hat[j*N+k]       = (r > (double) N/2) ? 0.0 : 1.0;
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
	int
	psi                  (nfft_plan& np) {
		
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
	int
	finalize             (nfft_plan& np, solver_plan_complex& spc) {
		
		solver_finalize_complex(&spc);
		nfft_finalize(&np);
		
		return 0;
		
	}
	
	
}

#endif //__NFFT_STUB_H__
