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
#endif
		
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
		extern int
		init                 (const int d, const int* N, const int M, const int* n, const int m, const nfft_plan* fnp, const solver_plan_complex* inp);
		
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
		extern int
		ift                  (const nfft_plan* np, const solver_plan_complex* spc, const int maxiter = 3, const double epsilon = 3e-7) ;
		
		/**
		 * @brief            Forward FT
		 *
		 * @param  np        NFFT plan
		 *
		 * @return           Success
		 */
		extern int
		ft                   (const nfft_plan* np);

		/**
		 * @brief            Adjoint FT
		 *
		 * @param  np        NFFT plan

		 * @return           Success
		 */
		extern int
		adjoint              (const nfft_plan* np);

		/**
		 * @brief            Set weights
		 * 
		 * @param  np        Plan
		 * @param  spc       Solver plan
		 * @return           Success
		 */
		extern int
		weights              (const nfft_plan* np, const solver_plan_complex* spc);

		/**
		 * @brief            Precompute PSI
		 *
		 * @param  np        NFFT plan
		 * @return           Success
		 */
		extern int
		psi                  (const nfft_plan* np);

		/**
		 * @brief            Finalise plans
		 *
		 * @param  np        Plan
		 * @param  spc       Solver plan
		 * @return           Success
		 */
		extern int
		finalize             (const nfft_plan* np, const solver_plan_complex* spc);

#ifdef __cplusplus
	}
#endif
	
}

#endif //__NFFT_STUB_H__
