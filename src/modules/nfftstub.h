#ifndef __NFFT_STUB_H__
#define __NFFT_STUB_H__

#include "nfft3util.h"
#include "nfft3.h"

/**
 * @brief Interface to <a href="http://www-user.tu-chemnitz.de/~potts/nfft/" target="NFFT">NFFT 3 library</a>
 */
namespace nfft {

#ifdef __cplusplus
	extern "C" {
#endif /* __cplusplus */

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
		 * @param  epsilon   Convergence criterium for solver (default: 3e-7)
		 *
		 * @return success
		 */
		extern int
		init                 (int d, int* N, int M, int* n, int m, nfft_plan* fnp, solver_plan_complex* inp, double epsilon  = 0.0000003);

		/**
		 * @brief            Inverse FT
		 * 
		 * @param  np        NFFT plan
		 * @param  spc       iNFFT plan
		 * @param  maxiter   Maximum NFFT iterations        
		 * @param  epsilon   Convergence criterium
		 *
		 * @return success
		 */
		extern int
		ift                  (nfft_plan* np, solver_plan_complex* spc, int maxiter = 3, double epsilon = 3e-7) ;
		
		/**
		 * @brief            Forward FT
		 *
		 * @param  np        NFFT plan
		 *
		 * @return success
		 */
		extern int
		ft                   (nfft_plan* np);

		/**
		 * @brief            Adjoint FT
		 *
		 * @param  np        NFFT plan
		 *
		 * @return success
		 */
		extern int
		adjoint              (nfft_plan* np);

		/**
		 * @brief            Set weights
		 */
		extern int
		weights              (nfft_plan* np, solver_plan_complex* spc);

		/**
		 * @brief            Precompute PSI
		 *
		 * @param  np        NFFT plan
		 */
		extern int
		psi                  (nfft_plan* np);

		/**
		 * @brief            Finalise
		 */
		extern int
		finalize             (nfft_plan* np, solver_plan_complex* spc);

#ifdef __cplusplus
	}
#endif /* __cplusplus */
	
}
#endif //__NFFT_STUB_H__
