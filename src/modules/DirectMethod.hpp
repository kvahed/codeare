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

#ifndef __DIRECT_METHOD_HPP__
#define __DIRECT_METHOD_HPP__

#include "ReconStrategy.hpp"

using namespace RRServer;


/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

	/**
	 * @brief Empty recon for test purposes
	 */
	class DirectMethod : public ReconStrategy {
		
		
	public:
		

		/**
		 * @brief Default constructor
		 */
		DirectMethod  () :
			m_np(1),
			m_verbose(true),
			m_dt(1.0) {};
		

		/**
		 * @brief Default destructor
		 */
		virtual 
		~DirectMethod () {

			this->Finalise();

		}
		
		
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
		
		double           m_dt;       /*!< @brief Simulation time steps                                        */
		bool             m_verbose;  /*!< @brief Verbose (Store magnetisation for every dt. HANDLE WITH CARE) */ 
		int              m_np;       /*!< @brief Number of OMP threads */
		int              m_ic;       /*!< @brief Perform intensity correction */

	};


	/**
	 * @brief         Time reversal for RF 
	 *
	 * @param  signal Acquired signal
	 * @param  pulse  Excitation pulse(s)
	 */
	void TimeReverseRF (const Matrix<cplx>& signal, const Matrix<double>& jac, Matrix<cplx>& pulse) {
	
		size_t nt = signal.Dim(0);
		size_t nc = signal.Dim(1);

		for (size_t i = 0; i < nt; i++)
			for (size_t c = 0; c < nc; c++)
				pulse(i,c) = (signal(nt-1-i,c)*(float)jac[nt-1-i]);
		
	}
	

	/**
	 * @brief         Time reversal for gradient trajectory
	 * 
	 * @param  acqgr  Acquisition gradients
	 * @param  excgr  Excitation gradients
	 */
	void TimeReverseGR (const Matrix<double>& acqgr, Matrix<double>& excgr) {

		size_t nt = acqgr.Dim(1);

		for (size_t i = 0; i < nt; i++) {
			
			excgr(0,i) = -acqgr(0,nt-1-i); 
			excgr(1,i) = -acqgr(1,nt-1-i); 
			excgr(2,i) = -acqgr(2,nt-1-i);

		}

	}


	/**
	 * @brief          Intensity correction
	 * 
	 * @param  b1maps  B1 maps
	 * @param  target  Target pattern
	 */
	void IntensityCorrection (const Matrix<cplx>& b1maps, Matrix<double>& target) {

		size_t nr = target.Dim(1);
		size_t nc = b1maps.Dim(1);
		float   a;

		cout << target.DimsToCString() << endl;
		cout << b1maps.DimsToCString() << endl;

		for (size_t i = 0; i < nr; i++) {

			a = 0.0;

			for (size_t j = 0; j < nc; j++)
				a += pow(abs(b1maps(i,j)),2);

			target(0,i) /= a;
			target(1,i) /= a;
			target(2,i) /= a;

		}

		target.MXDump("target.mat", "target");

	}

	/**
	 * @brief          CG NR algorithm
	 * 
	 * @param  maxit   Maximum # CG iterations
	 * @param  eps     Convergence criterium 
	 */
	/*
	void CGNR (// IN
			   const int& maxit, const float& eps, const Matrix<cplx>& rxm, const Matrix<cplx>& txm, Matrix<double> acqg, 

			   // OUT
			   Matrix<cplx>) {

		int           iters = 0;

		float         rn    = 0.0;
		float         an    = 0.0;

		raw           rtmp  = raw(0.0,0.0);
		
		vector<float> residue;

		// Convergence loop
		for (iters = 0; iters < maxit; iters++) {
			
			rn = pow(r.Norm().real(), 2);
			residue.push_back(rn/an);
			printf ("  %03i: CG residuum: %.9f\n", iters, residue.at(iters));
			
			// Convergence ? ----------------------------------------------
			if (std::isnan(residue.at(iters)) || residue.at(iters) <= eps)
				break;
			
			// Simulate Bloch receive mode
			Simulate (txm, rxm, rf, acqg,  r, target,  b0, m_dt, ACQUIRE, m_verbose, m_np, res, m);
			
			// Time reversal
			TimeReverseRF (res, j, rf);
			TimeReverseGR (ag, eg);

			// Simulate Bloch transmit mode
			Simulate (txp, rxm, rf, excg, sr, sample, sb0, m_dt,  EXCITE, m_verbose, m_np, res, m);
				
			rtmp  = (rn / (p.dotc(q)));
			a    += (p    * rtmp);
			s    += (stmp * rtmp);
			r    -= (q    * rtmp);
			p    *= cplx(pow(r.Norm().real(), 2)/rn);
			p    += r;
			
			// Verbose out put keeps all intermediate steps ---------------
			if (m_verbose) {
				image.Expand(m_dim); 
				signals.Expand(2);
				memcpy (  &image[(iters+1)*a.Size()], &a[0], a.Size()*sizeof(cplx));
				memcpy (&signals[(iters+1)*s.Size()], &s[0], s.Size()*sizeof(cplx));
			}
			
			
	}
	*/

}
#endif /* __DIRECT_METHOD_H__ */

