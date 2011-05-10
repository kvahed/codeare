#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include "nfft3util.h"
#include "nfft3.h"

int init   (int d, int* N, int M, int* n, int m, nfft_plan* fnp, solver_plan_complex* inp, double epsilon) {
	
	//assert (d == 2 || d == 3); 
	
	int      i, j, k, l;                            /* Some counters       */
	unsigned infft_flags = CGNR | PRECOMPUTE_DAMP;  /* flags for the infft */
	
#ifdef VERBOSE
	printf ("Initialising nfftplan with d:%i, Nx:%i, Ny:%i, M:%i, nx:%i, ny:%i\n", d, N[0], N[1], M, n[0], n[1]);
#endif
	
	nfft_init_guru (fnp, d, N, M, n, m,
					PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
					FFTW_MEASURE| FFTW_DESTROY_INPUT);
	
	if(fnp->nfft_flags & PRE_LIN_PSI)
		nfft_precompute_lin_psi(fnp);

	infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

	solver_init_advanced_complex (inp, (mv_plan_complex*)fnp, infft_flags);

	return 0;
	
}

int  weights (nfft_plan* np, solver_plan_complex* spc, double* w) {

	int j, k, z;

	memcpy (spc->w, w, np->M_total * sizeof(double));

	int N = np->N[0];

	if (spc->flags & PRECOMPUTE_DAMP) {
		if (np->d == 3) {
			for (j = 0; j < N; j++) 
				for (k = 0; k < N; k++)
					for (z = 0; z < N; z++) {
						int    j2 = j - N/2; 
						int    k2 = k - N/2;
						int    z2 = z - N/2;
						double r  = sqrt(j2*j2+k2*k2+z2*z2);
						spc->w_hat[z*N*N+j*N+k] = (r > (double) N/2) ? 0.0 : 1.0;
					}
		} else {
			for (j = 0; j < N; j++)
				for (k = 0; k < N; k++) {
					    int    j2 = j-N/2;
						int    k2 = k-N/2;
						double r  = sqrt(j2*j2+k2*k2);
						spc->w_hat[j*N+k]       = (r > (double) N/2) ? 0.0 : 1.0;
				}
		}
	}
	
	return np->M_total;

}


int  kspace  (nfft_plan* np, double* k) {

	memcpy (np->x, k, np->d * np->M_total * sizeof(double));
	return np->d * np->M_total;
	
}


int ft (nfft_plan* np, double* in, double* out) {

	memcpy (np->f_hat, in, 2 * np->N_total * sizeof (double));

	if(np->nfft_flags & PRE_PSI)
		nfft_precompute_psi (np);
	
	nfft_trafo (np);

	memcpy (out,    np->f, 2 * np->M_total * sizeof(double));

	return 0;

}


void ift (nfft_plan* np, solver_plan_complex* spc, double* in, double* out, int maxiter, double epsilon) {
	
	int j, k, l;

	memcpy (spc->y, in, 2 * np->M_total * sizeof(double));

	// precompute lin psi
	if(np->nfft_flags & PRE_PSI)
		nfft_precompute_psi(np);
	
	// precompute full psi 
	if(np->nfft_flags & PRE_FULL_PSI)
		nfft_precompute_full_psi(np);
	
	// init some guess 
	for(k = 0; k < np->N_total; k++)
		spc->f_hat_iter[k] = 0.0;
	
	double t = nfft_second();
	
	// inverse trafo 
	solver_before_loop_complex(spc);
   
	for( l = 0; l < maxiter; l++) {
		
		if (spc->dot_r_iter < epsilon) 
			break;

		solver_loop_one_step_complex(spc);

#ifdef VERBOSE		
		fprintf( stderr, "%e,  %i of %i\n", sqrt(spc->dot_r_iter), l+1, maxiter);
#endif		

	}

	t = nfft_second()-t;

#ifdef VERBOSE
	fprintf(stderr,"time: %e seconds.\n",t);
#endif

	memcpy (out, spc->f_hat_iter, 2 * np->N_total * sizeof (double));


}

int  finalize (nfft_plan* np, solver_plan_complex* spc) {

	solver_finalize_complex(spc);
	nfft_finalize(np);

}


