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
		
		omp_set_num_threads(NTHREADS);
		int tid      = omp_get_thread_num();
		
#pragma omp for
		for (int j = 0; j < ncoils; j++) {
			
			int    ipos   = j * imgsize;
			int    spos   = j * nsamples;
			raw    tmp    = raw(0.0, 0.0);

			// Copy data to FT
			for (int i = 0; i < imgsize; i++) {
				tmp = sm->at(ipos + i) * in->at(i);
				np[tid].f_hat[i][0] = tmp.real();
				np[tid].f_hat[i][1] = tmp.imag();
			}

			// Forward ft
			nfft::ft (&np[tid]);
			
			// Copy FTed data back
			for (int i = 0; i < nsamples; i++)
				out->at(i+spos) = raw(np[tid].f[i][0], np[tid].f[i][1]);
			
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
	int           ncoils   = sm->Dim(dim);
	int           nsamples = in->Size() / ncoils;
	int           imgsize  = out->Size();

	fftw_complex* ftout    = (fftw_complex*) malloc (imgsize * ncoils * sizeof (fftw_complex)); 

	// OMP Loop over coils, Inverse FT every signal in *in, 
	// Sum elementwise mutiplied images with according sensitivity maps 
#pragma omp parallel default (shared) 
	{
		
		omp_set_num_threads(NTHREADS);
		int tid      = omp_get_thread_num();
		
#pragma omp for

		for (int j = 0; j < ncoils; j++) {
			
			int    spos   = j * nsamples;
			int    ipos   = j * imgsize;
			
			// Copy to iFT
			for (int i = 0; i < nsamples; i++) {
				spc[tid].y[i][0] = (in->at(spos + i)).real();
				spc[tid].y[i][1] = (in->at(spos + i)).imag();
			}
			
			// Inverse FT
			nfft::ift (&np[tid], &spc[tid], maxit, epsilon);
			memcpy (&ftout[imgsize * j], &spc[tid].f_hat_iter[0][0], imgsize * sizeof (fftw_complex));

			
		}

		raw sens  = raw(0.0,0.0);
		int chunk = imgsize / NTHREADS;

#pragma omp for schedule (dynamic, chunk)

		for (int i = 0; i < imgsize; i++) 
			for (int j = 0; j < ncoils; j++) {
				int    ipos   = j * imgsize;
				sens        = sm->at(ipos + i);
				out->at(i) += raw(ftout[ipos + i][0], ftout[ipos + i][1]) * conj(sens);
			}

	}
	
	free (ftout);

	return OK;
	
}

