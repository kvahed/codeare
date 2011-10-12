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

#ifndef __COMPRESSED_SENSING_HPP__
#define __COMPRESSED_SENSING_HPP__

#include "ReconStrategy.hpp"
#include "FFT.hpp"
#include "DWT.hpp"

using namespace RRServer;


/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {


	struct CGParam {
		
		int    fft;
		int    pnorm;
		int    lsiter;
		int    lsto ;
		int    cgiter;
		
		double tvw;
		double xfmw;
		double cgconv;
		double l1smooth;
		double lsa;
		double lsb;
		
	};

	/**
	 * @brief CS reconstruction based on Sparse MRI v0.2 by Michael Lustig
	 */
	class CompressedSensing : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		CompressedSensing  () {};
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~CompressedSensing () {};
		
		
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
		Finalise () {
			return RRSModule::OK;
		};


	private:
		
		int            m_dim;    /**< Image recon dim */
		int            m_N[3];   /**< Data side lengths */
		int            m_csiter; /**< # global iterations */

		CGParam        m_cgparam;

	};


	void GradObj   (const Matrix<cplx>& x, const Matrix<cplx>& data, const CGParam& cgp, Matrix<cplx>& grad) {

		grad  = DWT::Backward (x);
		grad  = FFT::Forward  (x);
		grad -= data;
		grad  = FFT::Backward (x);
		grad  = DWT::Forward  (x);
		grad *= cplx(2.0,0.0);

	}

	void GradXFM   (const Matrix<cplx>& x, const Matrix<cplx>& data, const CGParam& cgp, Matrix<cplx>& grad) {
		
		Matrix<cplx> xtr = x.tr();

		grad  = x;
		grad *= xtr;
		grad += cplx(cgp.l1smooth,0.0);
		grad  = grad ^ (((float)cgp.pnorm)/2.0-1.0);
		grad *= x;

	}

	void GradTV    (const Matrix<cplx>& x, const Matrix<cplx>& data, const CGParam& cgp, Matrix<cplx>& grad) {

		// Dx = params.TV*(params.XFM'*x);
		// G = params.pNorm*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(params.pNorm/2-1);
		// grad = params.XFM*(params.TV'*G);

	}

	void WGradient          (const Matrix<cplx>& x, const Matrix<cplx>& data, const CGParam& cgp, Matrix<cplx>& grad) {

		GradObj (x, data, cgp, grad);
		GradXFM (x, data, cgp, grad);
		GradTV  (x, data, cgp, grad);

		//grad = (gradObj +  params.xfmWeight.*gradXFM + params.TVWeight.*gradTV);

	} 


	void NLCG               (Matrix<cplx>& x, Matrix<cplx>& data, const CGParam& cgp) {

		int k = 0;
		int t = 1;

		Matrix<cplx> grad (data.Dim());
		WGradient (x, data, cgp, grad);

	}


}
#endif /* __COMPRESSED_SENSING_H__ */
