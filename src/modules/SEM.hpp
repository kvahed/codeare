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

#include "OMP.hpp"

/**
 * @brief               Compute left hand side (i.e. Multiply E with spatial (image) data)
 *                      Forward NFFT in and elementwise multiply with spatial sensitivity of every channel 
 * 
 * @param  in           Original discretised sample O (Nx x Ny x Nz)
 * @param  sm           Sensitivity maps            O (Nx x Ny x Nz x Nc)
 * @param  np           Non-Cartesian strategy for non uniform ft
 * @param  out          Result                      O (Nk x Nc)
 * @param  dim          FT dimensions
 */
RRSModule::error_code 
E  (const Matrix<raw>& in, const Matrix<raw>& sm, nfft_plan* np, Matrix<raw>& out, const int& dim) {

	// Clear output container
	out.Zero();
	
	// Some dimensions
	int ncoils   = sm.Dim(dim);
	int nsamples = out.Size() / ncoils;
	int imgsize  = in.Size();

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
				tmp = sm[ipos + i] * in[i];
				np[tid].f_hat[i][0] = tmp.real();
				np[tid].f_hat[i][1] = tmp.imag();
			}

			// Forward ft
			nfft::ft (&np[tid]);
			
			// Copy FTed data back
			for (int i = 0; i < nsamples; i++)
				out[i+spos] = (raw(np[tid].f[i][0], np[tid].f[i][1]));
			
		}
		
	}
	
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
 * @param  dim          FT dimensions
 * @param  adjoint      Use adjoint?
 */
RRSModule::error_code
EH (const Matrix<raw>& in, const Matrix<raw>&  sm, nfft_plan*  np, solver_plan_complex* spc, const double& epsilon, 
	const int&      maxit,       Matrix<raw>& out, const int& dim, const bool&      adjoint) {

	// Clear outgoing container
	out.Zero();

	// Some dimensions
	int           ncoils   = sm.Dim(dim);
	int           nsamples = in.Size() / ncoils;
	int           imgsize  = out.Size();

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
				spc[tid].y[i][0] = (in[spos + i]).real();
				spc[tid].y[i][1] = (in[spos + i]).imag();
			}
			
			// Inverse FT or adjoint
			if (adjoint) {
				nfft::adjoint (&np[tid]);
				memcpy (&ftout[ipos], &np[tid].f_hat[0][0], imgsize*sizeof (fftw_complex));
			} else {
				nfft::ift (&np[tid], &spc[tid], maxit, epsilon);
				memcpy (&ftout[ipos], &spc[tid].f_hat_iter[0][0], imgsize*sizeof (fftw_complex));
			}
			
		}

		raw sens  = raw(0.0,0.0);
#pragma omp for schedule (dynamic, imgsize / NTHREADS)

		for (int i = 0; i < imgsize; i++) 
			for (int j = 0; j < ncoils; j++) {
				int pos = j * imgsize + i;
				sens = sm[pos];
				out[i] += (raw(ftout[pos][0], ftout[pos][1]) * conj(sens));
			}

	}
	
	free (ftout);

	return OK;
	
}

