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
#include "TVOP.hpp"

/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {


	/**
	 * @brief CG parameters
	 */
	struct CGParam {
		
		int    fft;
		int    pnorm;
		int    lsiter;
		int    lsto ;
		int    cgiter;
		
		double tvw;
		double xfmw;
		double cgconv;
		double l1;
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



	float Objective (const Matrix<cxfl>& ffdbx,       Matrix<cxfl>& ffdbg, const Matrix<cxfl>& ttdbx, 
					       Matrix<cxfl>& ttdbg, const Matrix<cxfl>&     x,       Matrix<cxfl>&     g, 
					 const Matrix<cxfl>&  data, const cxfl&             t,       float&         rmse,
					 const CGParam&        cgp) {
		
		float        obj = 0.0, tv = 0.0, xfm = 0.0;
		Matrix<cxfl> objm;
		float        p = (float)cgp.pnorm/2.0;

		objm  = ffdbx;
		objm += (ffdbg * t);
		objm -= data;
		obj   = pow(objm.Norm().real(), 2);

		
		if (cgp.tvw) {

			Matrix<cxfl> tvm;
			
			tvm  = ttdbx;
			tvm += (ttdbg * t);
			
			for (size_t i = 0; i < tvm.Size(); i++)
				tvm[i] *= conj(tvm[i]);
			
			tvm += cxfl(cgp.l1);
			tvm ^= p;
			
			tv   = 0.0;
			
			for (size_t i = 0; i < tvm.Size(); i++)
				tv += tvm[i].real();
			
		} 
		
		if (cgp.xfmw) {
			
			Matrix<cxfl> xfmm;
			
			xfmm  = x;
			xfmm += (g * t);
			
			for (size_t i = 0; i < xfmm.Size(); i++)
				xfmm[i] *= conj(xfmm[i]);
			
			xfmm += cxfl(cgp.l1);
			xfmm ^= p;
			
			xfm  = 0.0;
			
			for (size_t i = 0; i < xfmm.Size(); i++)
				xfm += xfmm[i].real();
			
		} 

		rmse = obj / pow(data.Norm().real(), 2);
		//printf (" RMS: %.9f\n", res.at(iters));

		return (obj + tv + xfm);

	}

	/**
	 * @brief Compute gradient of the data consistency
	 */
	Matrix<cxfl> GradObj (const Matrix<cxfl>& x, const Matrix<cxfl>& data, const CGParam& cgp) {
		
		Matrix<cxfl> g;
		
		g  = DWT::Backward (x);
		g  = FFT::Forward  (g);

		g -= data;

		g  = FFT::Backward (g);
		g  = DWT::Forward  (g);

		return g * cxfl(2.0,0.0);

	}


	/**
	 * @brief Compute gradient of L1-transform operator
	 */
	Matrix<cxfl> GradXFM   (const Matrix<cxfl>& x, const Matrix<cxfl>& data, const CGParam& cgp) {
		
		Matrix<cxfl> g = x;

		for (int i = 0; i < g.Size(); i++)
			g[i] *= conj(g[i]); 

		g += cxfl(cgp.l1);
		g ^= (((float)cgp.pnorm)/2.0-1.0);
		g *= x;

		return g * cxfl(cgp.xfmw);

	}


	/**
	 * @brief Compute gradient of the total variation operator
	 */
	Matrix<cxfl> GradTV    (const Matrix<cxfl>& x, const Matrix<cxfl>& data, const CGParam& cgp) {

		Matrix<cxfl> dx = TVOP::Transform(DWT::Backward(x));
		Matrix<cxfl> g = dx;

		for (int i = 0; i < g.Size(); i++)
			g[i] *= conj(g[i]); 

		g += cxfl(cgp.l1);
		g ^= (((float)cgp.pnorm)/2.0-1.0);

		for (int i = 0; i < g.Size(); i++)
			g[i] *= dx[i]; 

		g *= cxfl(cgp.pnorm);
		g  = TVOP::Adjoint (g);

		return g * cxfl(cgp.xfmw * cgp.tvw);

	}


	Matrix<cxfl> Gradient (const Matrix<cxfl>& x, const Matrix<cxfl>& data, const CGParam& cgp) {

		Matrix<cxfl> g = GradObj (x, data, cgp);
		
		if (cgp.xfmw)
			g += GradXFM (x, data, cgp);
		
		if (cgp.tvw)
			g += GradTV  (x, data, cgp);

		return g;

	} 


	Matrix<cxfl> NLCG     (const Matrix<cxfl>& x, const Matrix<cxfl>& data, const CGParam& cgp) {

		int   k  = 0;
		cxfl t0 = cxfl(1), t = cxfl(1);

		float rmse;

		Matrix<cxfl> g0, dx, ffdbx, ffdbg, ttdbx, ttdbg;

		g0 = Gradient (x, data, cgp);
		dx = -g0;
 
		ffdbx = FFT::Forward (DWT::Backward ( x));
		ffdbg = FFT::Forward (DWT::Backward (dx));
		
		if (cgp.tvw) {
			ttdbx = TVOP::Transform(DWT::Backward( x));
			ttdbg = TVOP::Transform(DWT::Backward(dx));
		}

		cxfl cz = cxfl(0.0);

		float f0 = Objective (ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data, cz, rmse, cgp);
		float f1 = Objective (ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data,  t, rmse, cgp);
		
		int i = 0;

		do {

			t *= cxfl(cgp.lsb);
			f1 = Objective(ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data, t, rmse, cgp);
			i++;

		} while (/*f1 > f0 - ((cxfl(cgp.lsa) * t) * Abs(g0.prodt(dx)) &&*/ i < cgp.lsiter);

			
		/*
lsiter = 0;

	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		[f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
	end

	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2
		t0 = t0 * beta;
	end 
	
	if lsiter<1
		t0 = t0 / beta;
	end

	x = (x + t*dx);

	%--------- uncomment for debug purposes ------------------------	
	disp(sprintf('%d   , obj: %f, RMS: %f, L-S: %d', k,f1,RMSerr,lsiter));

	%---------------------------------------------------------------
	
    %conjugate gradient calculation
    
	g1 = wGradient(x,params);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	%TODO: need to "think" of a "better" stopping criteria ;-)
	if (k > params.Itnlim) | (norm(dx(:)) < gradToll) 
		break;
	end
		*/

		return g0;

	}


}
#endif /* __COMPRESSED_SENSING_H__ */
