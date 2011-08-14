/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
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

#ifndef __RELATIVE_SENSITIVITIES_HPP__
#define __RELATIVE_SENSITIVITIES_HPP__

#include "ReconStrategy.hpp"

using namespace RRServer;


/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

	/**
	 * @brief Empty recon for test purposes
	 */
	class RelativeSensitivities : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		RelativeSensitivities  () {}
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~RelativeSensitivities () {}
		
		
		/**
		 * @brief Do nothing 
		 */
		virtual RRSModule::error_code
		Process ();
		
		/**
		 * @brief Do nothing 
		 */
		virtual RRSModule::error_code
		Init ();
		
		/**
		 * @brief Do nothing 
		 */
		virtual RRSModule::error_code
		Finalise ();



	private:

		
		double m_echo_shift;
		
		
	};

RRSModule::error_code
SVDCalibrate (const Matrix<raw>* imgs, Matrix<raw>* rxm, Matrix<raw>* txm, Matrix<double>* snro, Matrix<raw>* shim, const bool normalise) {

	int         nrxc = rxm->Dim(3);
	int         ntxc = txm->Dim(3);
	int      volsize = imgs->Dim(0) * imgs->Dim(1) * imgs->Dim(2);
	int       rtmsiz = nrxc * ntxc;
	int         vols = imgs->Size() / volsize / 2; // division by 2 (Echoes)
	int         rtms = imgs->Size() / rtmsiz / 2;  // division by 2 (Echoes)
	ticks        tic = getticks();
	
	printf ("  SVDing %i matrices of %lix%li ... ", rtms, nrxc, ntxc);
	
	// Permute dimensions on imgs for contiguous RAM access
	Matrix<raw> vxlm (nrxc, ntxc, imgs->Dim(0), imgs->Dim(1), imgs->Dim(2));
	
	for (int s = 0; s < imgs->Dim(2); s++)
		for (int l = 0; l < imgs->Dim(1); l++)
			for (int c = 0; c < imgs->Dim(0); c++)
				for (int t = 0; t < ntxc; t++)
					for (int r = 0; r < nrxc; r++)
						// multiplication with 2 (Need only 1st echo)
						vxlm (r, t, c, l, s) = imgs->At(c, l, s, 0, t, r); 
	
	Matrix<raw> OptSNR (imgs->Dim(0), imgs->Dim(1), imgs->Dim(2));
	int         threads = 1;
	
#pragma omp parallel default (shared) 
	{
		threads  = omp_get_num_threads();
	}
	
	Matrix<raw> m[threads]; // Combination
	Matrix<raw> u[threads]; // Left-side and  O(NRX x NRX)
	Matrix<raw> v[threads]; // Right-side singular vectors O(NTX, NTX)
	Matrix<raw> s[threads]; // Sorted singular values (i.e. first biggest) O (MIN(NRX,NTX));
	
	for (int i = 0; i < threads; i++) {
		m[i] = Matrix<raw>     (nrxc, ntxc);
		u[i] = Matrix<raw>     (nrxc, nrxc);
		v[i] = Matrix<raw>     (ntxc, ntxc);
		s[i] = Matrix<raw> (MIN(ntxc, nrxc), 1);
	}

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = rtms / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
		for (int i = 0; i < rtms; i++) {
			
			memcpy (&m[tid][0], &vxlm[i*rtmsiz], rtmsiz * sizeof(raw));
			
			m[tid].SVD ('A', &u[tid], &v[tid], &s[tid]);
			
			// U 
			for (int r = 0; r < nrxc; r++) rxm->At(r*volsize + i) = u[tid][r];

			// V is transposed!!! ------------------------------------------v
			for (int t = 0; t < nrxc; t++) txm->At(t*volsize + i) = v[tid][t*v[0].Dim(0)];
			  
			OptSNR[i] = real(s[tid][0]);
			
		}
		
	}

	printf ("done.                 (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
	
	return RRSModule::OK;

}


RRSModule::error_code
FTVolumes (Matrix<raw>* r) {
	
	long        imsize  = r->Dim(0) * r->Dim(1) * r->Dim(2);
	int         vols    = r->Size() / (imsize);
	int         threads = 1;
	Matrix<raw> hann    = Matrix<raw>::Ones(r->Dim(0), r->Dim(1), r->Dim(2)).HannWindow();
	ticks       tic     = getticks();
	
	printf ("  Fourier transforming %i volumes of %lix%lix%li ... ", vols, r->Dim(0), r->Dim(1), r->Dim(2));
	
#pragma omp parallel default (shared) 
	{
		threads  = omp_get_num_threads();
	}
	
	// # threads plans and matrices ------
	//
	// Note: FFTW plans should be created by 
	// one thread only!
	Matrix<raw> mr[threads];
	fftwf_plan 	 p[threads];
	
	for (int i = 0; i < threads; i++) {
		mr[i] = Matrix<raw>       (r->Dim(0), r->Dim(1), r->Dim(2));
		p[i]  = fftwf_plan_dft_3d (r->Dim(2), r->Dim(1), r->Dim(0), 
								   (fftwf_complex*)&mr[i][0], (fftwf_complex*)&mr[i][0], 
								   FFTW_BACKWARD, FFTW_ESTIMATE);
	}
	// ------------------------------------
	
	
	// ifftshift(ifftn(hann(fftshift(data))))
#pragma omp parallel default (shared)
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = vols / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
		for (int i = 0; i < vols; i++) {
			memcpy (&mr[tid][0], &r->At(i*imsize), imsize * sizeof(raw));
			mr[tid] = mr[tid].FFTShift();
			mr[tid] = mr[tid] * hann;
			fftwf_execute(p[tid]);
			mr[tid] = mr[tid].IFFTShift();
			memcpy (&r->At(i*imsize), &mr[tid][0], imsize * sizeof(raw));
		}

	}
	
	for (int i = 0; i < threads; i++)
		fftwf_destroy_plan(p[i]);

	printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

	
	return RRSModule::OK;
	
}


RRSModule::error_code 
RemoveOS (Matrix<raw>* imgs) {

	printf ("  Removing RO oversampling ... ");
	
	ticks tic    = getticks();
	Matrix<raw> tmp = (*imgs);
	
	imgs->Dim(0) = imgs->Dim(0) / 2;
	imgs->Reset();

	int   ossiz  = tmp.Dim(0);
	int   nssiz  = imgs->Dim(0);
	int   nscans = tmp.Size() / ossiz;
	int   offset = floor ( (float)(tmp.Dim(0)-imgs->Dim(0))/2.0 ); 

	for (int i = 0; i < nscans; i++)
		memcpy (&imgs->At(i*nssiz), &tmp.At(i*ossiz+offset), nssiz * sizeof(raw));

	printf ("done.                      (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
	
	return RRSModule::OK;

}


RRSModule::error_code
B0Map (const Matrix<raw>* imgs, Matrix<double>* b0, const float TE) {

	printf ("  Computing b0 maps ... ");

	ticks       tic = getticks();
	Matrix<raw> tmp;

	tmp = imgs->Mean(4);
	tmp.Squeeze();

	int nc = tmp.Dim(4);                           // Number of channels
	int np = tmp.Dim(0) * tmp.Dim(1) * tmp.Dim(2); // Number of pixels

	raw r;

#pragma omp parallel default (shared) private (r)
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = np / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
		for (int i = 0; i < np; i++) {
			r = raw (0.0,0.0);
			for (int j = 0; j < nc; j++)
				r += tmp[i + 2*j*np] * conj(tmp[i + (2*j+1)*np]);
			b0->At(i) = arg(r);// atan2 (r.imag(),r.real());
		}
	
	}

	printf ("done.                             (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

	return RRSModule::OK;

}




}
#endif /* __RELATIVE_SENSITIVITIES_H__ */

