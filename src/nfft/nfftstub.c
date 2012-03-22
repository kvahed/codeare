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

#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include "nfft3util.h"
#include "nfft3.h"

int init   (int d, int* N, int M, int* n, int m, nfft_plan* np, solver_plan_complex* inp) {
	
	unsigned infft_flags = CGNR | PRECOMPUTE_DAMP;  /* flags for the infft */

#ifdef VERBOSE
	printf ("Initialising nfftplan with d:%i, Nx:%i, Ny:%i, Nz:%i, M:%i, nx:%i, ny:%i, nz:%i\n", d, N[0], N[1], N[2], M, n[0], n[1], n[2]);
#endif
	
	nfft_init_guru (np, d, N, M, n, m,
					PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
					FFTW_MEASURE| FFTW_DESTROY_INPUT);
	
	if(np->nfft_flags & PRE_LIN_PSI)
		nfft_precompute_one_psi(np);

	infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

	solver_init_advanced_complex (inp, (mv_plan_complex*)np, infft_flags);

	return 0;
	
}


int weights (nfft_plan* np, solver_plan_complex* spc) {
	
	int j, k, z, N = np->N[0];
	
	if (spc->flags & PRECOMPUTE_DAMP) {
		if (np->d == 3) {
			for (j = 0; j < N; j++) {
				int    j2 = j - N/2; 
				for (k = 0; k < N; k++) {
					int    k2 = k - N/2;
					for (z = 0; z < N; z++) {
						int    z2 = z - N/2;
						double r  = sqrt(j2*j2+k2*k2+z2*z2);
						spc->w_hat[z*N*N+j*N+k] = (r > (double) N/2) ? 0.0 : 1.0;
					}
				}
			}
		} else {
			for (j = 0; j < N; j++) {
				int    j2 = j-N/2;
				for (k = 0; k < N; k++) {
					int    k2 = k-N/2;
					double r  = sqrt(j2*j2+k2*k2);
					spc->w_hat[j*N+k]       = (r > (double) N/2) ? 0.0 : 1.0;
				}
			}
		}
	}

	return np->M_total;
	
}
	

int ft (nfft_plan* np) {

	nfft_trafo (np);
	return 0;

}


int adjoint (nfft_plan* np) {

	nfft_adjoint (np);
	return 0;

}


void psi (nfft_plan* np) {

	/* precompute full psi */
	if(np->nfft_flags & PRE_PSI)
		nfft_precompute_one_psi(np);

	/* precompute full psi */
	if(np->nfft_flags & PRE_FULL_PSI)
		nfft_precompute_full_psi(np);
	
}

void ift (nfft_plan* np, solver_plan_complex* spc, int maxiter, double epsilon) {
	
	int k, l, err = 0;
	
	/* init some guess */
	for(k = 0; k < np->N_total; k++)
		spc->f_hat_iter[k] = 0.0;
	
	double t = nfft_second();
	
	/* inverse trafo */
	solver_before_loop_complex(spc);
	
	for (l = 0; l < maxiter; l++) {
		
		if (spc->dot_r_iter < epsilon) 
			break;
		
		solver_loop_one_step_complex(spc);
		
#ifdef NVERBOSE		
		fprintf( stderr, "%e,  %i of %i\n", sqrt(spc->dot_r_iter), l+1, maxiter);
#endif		
		
	}
	
	t = nfft_second()-t;
	
#ifdef NVERBOSE
	fprintf(stderr,"nfft time: %.4f seconds.\n",t);
#endif
	
}

int  finalize (nfft_plan* np, solver_plan_complex* spc) {

	solver_finalize_complex(spc);
	nfft_finalize(np);

	return 0;

}


