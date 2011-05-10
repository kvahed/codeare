#ifndef __NFFT_STUB_H__
#define __NFFT_STUB_H__

#define GAMMA 4.15

#include "nfft3util.h"
#include "nfft3.h"

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
		 * @param  epsilon   Convergence criterium for solver (default: 3e-7)
		 * @parma  fnp       Forward FT plan
		 * @param  inp       Inverse FT plan
		 * @param  scp       Solver plan
		 */
		extern const int
		init                 (int d, int* N, int M, int* n, int m, nfft_plan* fnp, solver_plan_complex* inp, double epsilon  = 0.0000003);

		/**
		 * @brief            Inverse FT
		 * 
		 * @param  in        Data for inverse FT
		 * @param  k         K-space trajectory
		 * @param  dk        Jacobian determinants of k(t). (i.e. det(D_f(k))
		 * @param  out       FTed data
		 */
		extern const int
		ift                  (nfft_plan* np, solver_plan_complex* spc, double* in, double* out, int maxiter = 3, double epsilon = 3e-7) ;
		
		/**
		 * @brief            Forward FT
		 *
		 * @param  in        Data for forward FT
		 * @param  k         K-space trajectory
		 * @param  out       FTed data 
		 */
		extern int
		ft                   (nfft_plan* np, double* in, double* out);

		extern const int
		weights              (nfft_plan* np, solver_plan_complex* spc, double* w);

		extern const int
		kspace               (nfft_plan* np, double* k);

		extern const int
		finalize             (nfft_plan* np, solver_plan_complex* spc);

#ifdef __cplusplus
	}
#endif /* __cplusplus */
	
}
#endif //__NFFT_STUB_H__
