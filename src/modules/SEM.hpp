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
#include "Access.hpp"

template <class T> inline static Matrix<T> 
E (const Matrix<T>& in, const Matrix<T>& sm, NFFT<T>** fts, const int& dim) {

	// Some dimensions
	size_t nc, nr, nk;

	nc = size  (sm,dim); 
	nr = numel (in);
	nk = fts[0]->FPlan()->M_total;

	Matrix<T> res (nk,nc);
	
	// Loop over coils, Elementwise multiplication of maps with in (s.*in), ft and store in out
	
#pragma omp parallel default (shared) 
	{

		NFFT<T>&  ft = *fts[omp_get_thread_num()];
		Matrix<T> tmp;
		
#pragma omp for // coils
		for (int j = 0; j < nc; j++) {
			
			tmp  = (dim == 2) ? Slice (sm, j) : Volume (sm, j);
			tmp *= in;
			tmp  = ft * tmp;
			
			memcpy (&res[j*nk], &tmp[0], nk * sizeof(T));
			
		}
		
	}
	
	return res;
	
}


template <class T> inline static Matrix<T> 
EH (const Matrix<T>& in, const Matrix<T>& sm, NFFT<T>** fts, const int& dim, 
	const double& epsilon = 7.0e-3, const int& maxit = 1, const bool& adjoint = false) {

	// Some dimensions
	size_t nr, nc, nk;

	nc = size(sm,dim);
	nk = size(in,0);
	nr = numel(sm)/nc;

	Matrix<T> res 
		(size(sm,0), 
		 size(sm,1),
		 (dim == 3) ? size(sm,2) : nc, 
		 (dim == 3) ?         nc :  1);

	// OMP Loop over coils, Inverse FT every signal in *in, 
	// Sum elementwise mutiplied images with according sensitivity maps 
#pragma omp parallel default (shared) 
	{
		
		NFFT<T>&  ft = *fts[omp_get_thread_num()];
		Matrix<T> tmp;
		
#pragma omp for
		for (int j = 0; j < nc; j++) {
			
			int    spos   = j * nk;
			int    ipos   = j * nr;
			
			if      (dim == 2)  Slice (res, j, ft ->* Column (in, j) *  Slice (sm, j));
			else if (dim == 3) Volume (res, j, ft ->* Column (in, j) * Volume (sm, j));
			
		}
		
	}
	
	return sum (res, dim);

}




/**
 * @brief               Compute left hand side (i.e. Multiply E with spatial (image) data)
 *                      Forward NFFT in and elementwise multiply with spatial sensitivity of every channel 
 * 
 * @param  in           Original discretised sample O (Nx x Ny x Nz)
 * @param  sm           Sensitivity maps            O (Nx x Ny x Nz x Nc)
 * @param  ic           Intesity correction map
 * @param  np           Non-Cartesian strategy for non uniform ft
 * @param  out          Result                      O (Nk x Nc)
 * @param  dim          FT dimensions
 */
/*
static RRSModule::error_code 
E  (const Matrix<cxdb>& in, const Matrix<cxdb>& sm, const Matrix<double>& ic, NFFT<cxdb>** np, Matrix<cxdb>& out, const int& dim) {

	// Some dimensions
	int ncoils   = size(sm, dim);
	int nsamples = numel(out) / ncoils;
	int imgsize  = numel(in);


	// Loop over coils, Elementwise multiplication of maps with in (s.*in), ft and store in out
	
#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		
#pragma omp for 
		for (int i = 0; i < numel(ic); i++)
			in[i] *= ic[i]; 
		
#pragma omp for
		for (int j = 0; j < ncoils; j++) {
			
			int    ipos   = j * imgsize;
			int    spos   = j * nsamples;
			raw    tmp    = cxdb (0.0, 0.0);

			// Copy data to FT
			for (int i = 0; i < imgsize; i++) {

				tmp = sm[ipos + i] * in[i];

				np[tid]->FPlan->f_hat[i][0] = tmp.real();
				np[tid]->FPlan->f_hat[i][1] = tmp.imag();

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

*/



/**
 * @brief               Compute left hand side (i.e. Multiply E with spatial (image) data)
 *                      Forward NFFT in and elementwise multiply with spatial sensitivity of every channel 
 * 
 * @param  in           Original discretised sample O (Nx x Ny x Nz)
 * @param  sm           Sensitivity maps            O (Nx x Ny x Nz x Nc)
 * @param  ic           Intesity correction map
 * @param  np           Non-Cartesian strategy for non uniform ft
 * @param  out          Result                      O (Nk x Nc)
 * @param  dim          FT dimensions
 */
RRSModule::error_code 
E  (const Matrix<cxfl>& in, const Matrix<cxfl>& sm, const Matrix<double>& ic, nfft_plan* np, Matrix<cxfl>& out, const int& dim) {

	// Some dimensions
	int ncoils   = sm.Dim(dim);
	int nsamples = out.Size() / ncoils;
	int imgsize  = in.Size();


	// Loop over coils, Elementwise multiplication of maps with in (s.*in), ft and store in out
	
#pragma omp parallel default (shared) 
	{

		omp_set_num_threads(8);
		int tid      = omp_get_thread_num();

#pragma omp for 
		for (int i = 0; i < ic.Size(); i++)
			in[i] *= ic[i]; 
		
#pragma omp for
		for (int j = 0; j < ncoils; j++) {
			
			int    ipos   = j * imgsize;
			int    spos   = j * nsamples;
			cxfl    tmp    = cxfl(0.0, 0.0);

			// Copy data to FT
			for (int i = 0; i < imgsize; i++) {

				tmp = sm[ipos + i] * in[i] ;

				np[tid].f_hat[i][0] = tmp.real();
				np[tid].f_hat[i][1] = tmp.imag();

			}

			// Forward ft
			nfft::ft (&np[tid]);
			
			// Copy FTed data back
			for (int i = 0; i < nsamples; i++)
				out[i+spos] = (cxfl(np[tid].f[i][0], np[tid].f[i][1]));
			
		}
		
	}
	
	return OK;
	
}

/**
 * @brief               Compute right hand side (i.e. Multiply E^H, Hermitian counterpart to E, with k-space data)
 *
 * @param  in           K-space samples along trajectory O (Nk x Nc)
 * @param  sm           Sensitivity maps                 O (Nx x Ny x Nz x Nc)
 * @param  ic           Intesity correction map
 * @param  np           NuFFT plan
 * @param  spc          Solver plan
 * @param  epsilon      Convergence criterium for ift (default 3e-7)
 * @param  maxit        Maximum number of solver iterations (default 3)
 * @param  out          Returned product                 O (Nx x Ny x Nz)
 * @param  dim          FT dimensions
 * @param  adjoint      Use adjoint?
 */
RRSModule::error_code
EH (const Matrix<cxfl>& in, const Matrix<cxfl>&  sm, const Matrix<double>& ic, 
	nfft_plan* np, solver_plan_complex* spc, const double& epsilon, const int& maxit,
	Matrix<cxfl>& out, const int& dim, const bool& adjoint) {

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
		
		omp_set_num_threads(8);
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
				memcpy (&ftout[ipos],       &np[tid].f_hat[0][0], imgsize*sizeof (fftw_complex));
			} else {
				nfft::ift (&np[tid],  &spc[tid], maxit, epsilon);
				memcpy (&ftout[ipos], &spc[tid].f_hat_iter[0][0], imgsize*sizeof (fftw_complex));
			}
			
		}

		cxfl sens  = cxfl(0.0,0.0);

#pragma omp for schedule (guided, 1024)
		for (int i = 0; i < imgsize; i++) {

			for (int j = 0; j < ncoils; j++) {
				int pos = j * imgsize + i;
				out[i] += (cxfl(ftout[pos][0], ftout[pos][1]) * conj(sm[pos]));
			}

			out[i] *= ic[i];

		}
	}
	

	free (ftout);

	return OK;
	
}

