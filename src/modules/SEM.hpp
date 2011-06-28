#include "OMP.hpp"

/**
 * @brief               Compute left hand side (i.e. Multiply E with spatial (image) data)
 *                      Forward NFFT in and elementwise multiply with spatial sensitivity of every channel 
 * 
 * @param  in           Original discretised sample O (Nx x Ny x Nz)
 * @param  sm           Sensitivity maps            O (Nx x Ny x Nz x Nc)
 * @param  np           Non-Cartesian strategy for non uniform ft
 * @param  out          Result                      O (Nk x Nc)
 */
RRSModule::error_code 
E  (Matrix<raw>* in, Matrix<raw>* sm, nfft_plan* np, Matrix<raw>* out, int dim) {

	// Clear output container
	out->Zero();
	
	// Some dimensions
	int        ncoils   = sm->Dim(dim);
	int        nsamples = out->Size() / ncoils;
	int        imgsize  = in->Size();

	// Loop over coils, Elementwise multiplication of maps with in (s.*in), ft and store in out
	
#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		omp_set_num_threads(NTHREADS);
		
#pragma omp for
		for (int j = 0; j < ncoils; j++) {
			
			int    ipos   = j * imgsize;
			int    spos   = j * nsamples;
			double dt [2] = {0.0, 0.0};
			raw    tmp    = raw(0.0, 0.0);
			int    sof    = sizeof (fftw_complex);

			// Copy data to FT
			for (int i = 0; i < imgsize; i++) {
				
				tmp = sm->at(ipos + i) * in->at(i);
				dt[0] = tmp.real();
				dt[1] = tmp.imag();
				memcpy ( &(np[tid].f_hat)[i], dt, sof);
				
			}
			
			// Forward ft
			nfft::ft (&np[tid]);
			
			// Copy FTed data back
			for (int i = 0; i < nsamples; i++) {
				memcpy (dt, &(np[tid].f)[i+spos], sof);
				out->at(i+spos) = raw(dt[0],dt[1]);
			}
			
		}
		
	}
	
	// Return success
	return OK;
	
}

/**
 * @brief               Compute right hand side (i.e. Multiply E^H, Hermitian counterpart to E, with k-space data)
 *
 * @param  in           K-space samples along trajectory O (Nk x Nc)
 * @param  sm           Sensitivity maps                 O (Nx x Ny x Nz x Nc)
 * @param  np           NuFFT plan
 * @param  spc          Solver plan
 * @param  epsilon      Convergence criterium for ift (default 3e-7)
 * @param  maxit        Maximum number of solver iterations (default 3)
 * @param  out          Returned product                 O (Nx x Ny x Nz)
 */
RRSModule::error_code
EH (Matrix<raw>* in, Matrix<raw>* sm, nfft_plan* np, solver_plan_complex* spc, double epsilon, int maxit, Matrix<raw>* out, int dim) {

	// Clear outgoing container
	out->Zero();

	// Some dimensions
	int        ncoils   = sm->Dim(dim);
	int        nsamples = in->Size() / ncoils;
	int        imgsize  = out->Size();

	double* ftout = (double*) malloc (2 * ncoils * imgsize  * sizeof(double));
	
	// OMP Loop over coils, Inverse FT every signal in *in, 
	// Sum elementwise mutiplied images with according sensitivity maps 
#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		omp_set_num_threads(NTHREADS);
		
#pragma omp for
		for (int j = 0; j < ncoils; j++) {
			
			// Containers for FT I/O
			double* ftin  = (double*) malloc (2 * nsamples * sizeof(double));
			
			int    spos   = j * nsamples;
			int    ipos   = j * imgsize;
			
			// Copy to iFT
			for (int i = 0; i < nsamples; i++) {
				ftin[2*i  ] = (in->at(spos + i)).real();
				ftin[2*i+1] = (in->at(spos + i)).imag();
			}
			
			// Inverse FT
			nfft::ift (&np[tid], &spc[tid], ftin, &ftout[2*ipos], maxit, epsilon);
			
			//free (ftin);
			
		}
	}

	int chunk    = out->Size()/NTHREADS; 
		
#pragma omp for schedule(dynamic,chunk)

	for (int i = 0; i < out->Size(); i++)
		for (int j = 0; j < ncoils; j++) {
			int ipos = j * imgsize;
			raw sens = sm->at(ipos + i);
			out->at(i) += raw(ftout[2*i+2*ipos], ftout[2*i+1+2*ipos]) * conj(sens);
		}
	
	// Free RAM
	free (ftout);

	return OK;
	
}

