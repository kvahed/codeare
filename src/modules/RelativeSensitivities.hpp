/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Daniel Brenner
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

const static float GAMMA_1_PER_UT_MS = 2.675222099;

	/**
	 * @brief b0 abd Relative b1 maps from 2SPGREs 
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
		double m_cutoff;

		int m_use_bet;
		int m_log_mask;
		int m_weigh_maps;
		
		
	};

	RRSModule::error_code
	SVDCalibrate (const Matrix<cplx>& imgs, Matrix<cplx>& rxm, Matrix<cplx>& txm, Matrix<double>& snro, Matrix<cplx>& shim, const bool& normalise) {
		
		size_t    nrxc = rxm.Dim(3);
		size_t    ntxc = txm.Dim(3);
		size_t volsize = imgs.Dim(0) * imgs.Dim(1) * imgs.Dim(2);
		size_t  rtmsiz = nrxc * ntxc;
		size_t    vols = imgs.Size() / volsize / 2; // division by 2 (Echoes)
		size_t    rtms = imgs.Size() / rtmsiz / 2;  // division by 2 (Echoes)
		ticks      tic = getticks();

		printf ("  SVDing %i matrices of %ix%i ... ", (int)rtms, (int)nrxc, (int)ntxc); fflush(stdout);
		
		// Permute dimensions on imgs for contiguous RAM access
		Matrix<cplx> vxlm (nrxc, ntxc, imgs.Dim(0), imgs.Dim(1), imgs.Dim(2));
		
		for (size_t s = 0; s < imgs.Dim(2); s++)
			for (size_t l = 0; l < imgs.Dim(1); l++)
				for (size_t c = 0; c < imgs.Dim(0); c++)
					for (size_t t = 0; t < ntxc; t++)
						for (size_t r = 0; r < nrxc; r++)
							// multiplication with 2 (Need only 1st echo)
							vxlm(r, t, c, l, s) = imgs(c, l, s, 0, t, r); 
		
		int         threads = 1;
		
#pragma omp parallel default (shared) 
		{
			threads  = omp_get_num_threads();
		}
		
		Matrix<cplx> m[threads]; // Combination
		Matrix<cplx> u[threads]; // Left-side and  O(NRX x NRX)
		Matrix<cplx> v[threads]; // Right-side singular vectors O(NTX, NTX)
		Matrix<cplx> s[threads]; // Sorted singular values (i.e. first biggest) O (MIN(NRX,NTX));
		
		for (int i = 0; i < threads; i++) {
			m[i] = Matrix<cplx>     (nrxc, ntxc);
			u[i] = Matrix<cplx>     (nrxc, nrxc);
			v[i] = Matrix<cplx>     (ntxc, ntxc);
			s[i] = Matrix<cplx> (MIN(ntxc, nrxc), 1);
		}

#pragma omp parallel default (shared) 
		{
			
			int tid = omp_get_thread_num();
			
#pragma omp for schedule (dynamic, rtms / omp_get_num_threads())
			
			for (size_t i = 0; i < rtms; i++) {
				
				memcpy (&m[tid][0], &vxlm[i*rtmsiz], rtmsiz * sizeof(cplx));
				
				m[tid].SVD (&u[tid], &v[tid], &s[tid], 'S');
				
				for (size_t r = 0; r < nrxc; r++) rxm[r*volsize + i] = u[tid][r]             * exp(cplx(0.0,1.0)*arg(u[tid][0]));             // U 
				for (size_t t = 0; t < nrxc; t++) txm[t*volsize + i] = v[tid][t*v[0].Dim(0)] * exp(cplx(0.0,1.0)*arg(v[tid][0])); // V is transposed!!!
				
				snro[i] = real(s[tid][0]);
				
			}
			
		}
		
		printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

		return RRSModule::OK;
		
	}
	
	
	RRSModule::error_code
	FTVolumes (Matrix<cplx>& r) {
		
		long         imsize  = r.Dim(0) * r.Dim(1) * r.Dim(2);
		size_t       vols    = r.Size() / (imsize);
		int          threads = 1;
		
		Matrix<cplx> hann (r.Dim(0), r.Dim(1), r.Dim(2));
		hann = cplx(1.0,0.0);
		hann = hann.HannWindow();
						   
		ticks        tic     = getticks();
		
		printf ("  Fourier transforming %i volumes of %ix%ix%i ... ", (int)vols, (int)r.Dim(0), (int)r.Dim(1), (int)r.Dim(2)); fflush(stdout);

#pragma omp parallel default (shared) 
		{
			threads  = omp_get_num_threads();
		}
		
		// # threads plans and matrices ------
		//
		// Note: FFTW plans should be created by 
		// one thread only!
		Matrix<cplx> mr[threads];
		fftwf_plan 	  p[threads];
		
		for (int i = 0; i < threads; i++) {
			mr[i] = Matrix<cplx>      (r.Dim(0), r.Dim(1), r.Dim(2));
			p[i]  = fftwf_plan_dft_3d (r.Dim(2), r.Dim(1), r.Dim(0), 
									   (fftwf_complex*)&mr[i][0], (fftwf_complex*)&mr[i][0], 
									   FFTW_BACKWARD, FFTW_ESTIMATE);
		}
		// ------------------------------------
		
		
#pragma omp parallel default (shared)
		{
			
			int tid = omp_get_thread_num();
			
#pragma omp for schedule (dynamic, vols / omp_get_num_threads())
			
			for (size_t i = 0; i < vols; i++) {
				memcpy (&mr[tid][0], &r[i*imsize], imsize * sizeof(cplx));
				mr[tid]  = mr[tid].FFTShift();
				mr[tid]  = mr[tid] * hann;
				fftwf_execute(p[tid]);
				mr[tid] = mr[tid].IFFTShift();
				memcpy (&r[i*imsize], &mr[tid][0], imsize * sizeof(cplx));
			}
			
		}
		
		for (int i = 0; i < threads; i++)
			fftwf_destroy_plan(p[i]);
		
		printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
		
		return RRSModule::OK;
		
	}
	
	RRSModule::error_code 
	RemoveOS (Matrix<cplx>& imgs) {
		
		printf ("  Removing RO oversampling ... "); fflush(stdout);
	
		ticks tic = getticks();
		Matrix<cplx> tmp = imgs;
		
		imgs.Dim(0) = imgs.Dim(0) / 2;
		imgs.Reset();
		
		size_t   ossiz  = tmp.Dim(0);
		size_t   nssiz  = imgs.Dim(0);
		size_t   nscans = tmp.Size() / ossiz;
		size_t   offset = floor ( (float)(tmp.Dim(0)-imgs.Dim(0))/2.0 ); 
		
		for (size_t i = 0; i < nscans; i++)
			memcpy (&imgs[i*nssiz], &tmp[i*ossiz+offset], nssiz * sizeof(cplx));
		
		printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
		
		return RRSModule::OK;
		
	}


	RRSModule::error_code
	B0Map (const Matrix<cplx>& imgs, Matrix<double>& b0, const float& TE) {
		
		printf ("  Computing b0 maps ... "); fflush(stdout);
		
		ticks  tic = getticks();
		Matrix<cplx> tmp;
		
		tmp = imgs.Mean(4);
		tmp.Squeeze();
		
		size_t nc = tmp.Dim(4);                           // Number of channels
		size_t np = tmp.Dim(0) * tmp.Dim(1) * tmp.Dim(2); // Number of pixels
		
		cplx   r;
		
#pragma omp parallel default (shared) private (r)
		{
			
#pragma omp for schedule (dynamic, np / omp_get_num_threads())
			
			for (int i = 0; i < np; i++) {
				r = cplx (0.0,0.0);
				for (int j = 0; j < nc; j++)
					r += tmp[i + 2*j*np] * conj(tmp[i + (2*j+1)*np]);
				b0[i] = arg(r) / GAMMA_1_PER_UT_MS / TE / 2; 

			}
			
		}
		
		printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
		
		return RRSModule::OK;
		
	}


	/**
	 * @brief Brain segmentation with FSL bet
	 *
	 * @brief In-out Images/mask
	 */
	int SegmentBrain (const Matrix<double>& img, Matrix<short>& msk) {
		
		printf ("  Brain segmentation with FSL(bet2) ... "); fflush(stdout);

		ticks  tic = getticks();

		std::string orig = "orig.nii.gz";
		std::string mask = "mask_mask.nii.gz";
		std::string cmd  = "/usr/local/bin/mask.sh";

		img.NIDump(orig);
		printf ("exporting ... "); fflush(stdout);
		
		if (!std::system(NULL))
			return FATAL_SYSTEM_CALL;

		// Call bet
		printf ("bet2ing ... "); fflush(stdout);
		system (cmd.c_str());

		printf ("importing ... "); fflush(stdout);
		msk.NIRead(mask);

 		printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

		return RRSModule::OK;

	}



	void LogMask (Matrix<double>& m, const double& m_cutoff) {

		for (int i = 0; i < m.Size(); i++)
			m[i] = (log(abs(m[i])) < m_cutoff) ? 0.0 : 1.0;

	}


}
#endif /* __RELATIVE_SENSITIVITIES_H__ */

