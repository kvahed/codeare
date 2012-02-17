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
#include "CX.hpp"

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
					       Matrix<cxfl>& ttdbg,       Matrix<cxfl>&     x,       Matrix<cxfl>&     g, 
					 const Matrix<cxfl>&  data, const cxfl&             t,       float&         rmse,
					 const CGParam&        cgp) {
		
		float        obj = 0.0, tv = 0.0, xfm = 0.0;
		Matrix<cxfl> objm;
		float        p = (float)cgp.pnorm/2.0;

		objm  =  ffdbx;
		if (cabs(t) > 0.0)
			objm += (ffdbg * t);
		objm -=  data;
		obj   =  creal(objm.dotc(objm));

		if (cgp.tvw) {

			Matrix<cxfl> tvm;
			
			tvm  =  ttdbx;
			if (cabs(t) > 0.0)
				tvm += (ttdbg * t);
			tvm *=  CX::Conj(tvm);
			tvm +=  cxfl(cgp.l1);
			tvm ^=  p;
			
			for (size_t i = 0; i < tvm.Size(); i++)
				tv += tvm[i].real();
			
		} 

		if (cgp.xfmw) {
			
			Matrix<cxfl> xfmm;
			
			xfmm  =  x;
			if (cabs(t) > 0.0)
				xfmm += (g * t);
			xfmm *=  CX::Conj(xfmm);
			xfmm +=  cxfl(cgp.l1);
			xfmm ^=  p;
			
			for (size_t i = 0; i < xfmm.Size(); i++)
				xfm += xfmm[i].real();
			
		} 

		float sda = 0.0;
		
		for (size_t i = 0; i < data.Size(); i++)
			if (cabs(data[i])>0)
				sda += cabs(data[i]);
		
		rmse = sqrt(obj / sda);

		return (obj + tv + xfm);

	}

	/**
	 * @brief Compute gradient of the data consistency
	 */
	Matrix<cxfl> GradObj ( Matrix<cxfl>& x, const Matrix<cxfl>& data, const CGParam& cgp) {
		
		Matrix<cxfl> g = x;
		
		g  = FFT::Forward (DWT::Backward (g));
		g -= data;
		g  = DWT::Forward (FFT::Backward (g));
		g *= cxfl(2.0);

		g.MXDump ("gfun.mat", "g");

		return g;

	}


	/**
	 * @brief Compute gradient of L1-transform operator
	 */
	Matrix<cxfl> GradXFM   (Matrix<cxfl>& x, const CGParam& cgp) {
		
		Matrix<cxfl> g;

		g  = CX::Conj(x) * x;
		g += cxfl(cgp.l1);
		g ^= (((float)cgp.pnorm)/2.0-1.0);
		g *= x * cgp.xfmw;

		g.MXDump ("gxfm.mat", "gxfm");

		return g;

	}


	/**
	 * @brief Compute gradient of the total variation operator
	 */
	Matrix<cxfl> GradTV    (const Matrix<cxfl>& x, const CGParam& cgp) {

		Matrix<cxfl> dx = TVOP::Transform(DWT::Backward(x));
		Matrix<cxfl> g  = CX::Conj(dx) * dx ;

		g += cxfl(cgp.l1);
		g ^= (((float)cgp.pnorm)/2.0-1.0);

		for (int i = 0; i < g.Size(); i++)
			g[i] *= dx[i]; 

		g *= cxfl(cgp.pnorm);
		g  = DWT::Forward(TVOP::Adjoint(g));
		g *= cxfl(cgp.tvw);

		g.MXDump ("gtv.mat", "gtv");

		return g;

	}


	Matrix<cxfl> Gradient (Matrix<cxfl>& x, const Matrix<cxfl>& data, const CGParam& cgp) {

		Matrix<cxfl> g = GradObj (x, data, cgp);
		
		if (cgp.xfmw)
			g += GradXFM (x, cgp);
		
		if (cgp.tvw)
			g += GradTV  (x, cgp);

		return g;

	} 


	void NLCG     (Matrix<cxfl>& x, const Matrix<cxfl>& data, const CGParam& cgp) {

		int          k  = 0;
		cxfl         t0 = cxfl(1), t = cxfl(1);
		Matrix<cxfl> x0 = x;

		float rmse, bk;

		Matrix<cxfl> g0, g1, dx, ffdbx, ffdbg, ttdbx, ttdbg;

		g0 = Gradient (x, data, cgp);
		dx = -g0;

		do {

			t = t0;

			ffdbx = FFT::Forward (DWT::Backward ( x));
			ffdbg = FFT::Forward (DWT::Backward (dx));
			
			if (cgp.tvw) {
				ttdbx = TVOP::Transform (DWT::Backward( x));
				ttdbg = TVOP::Transform (DWT::Backward(dx));
			}
			
			cxfl cz = cxfl(0.0);
			
			float f0 = Objective (ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data, cz, rmse, cgp);
			float f1 = Objective (ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data,  t, rmse, cgp);
			
			int i = 0;
			
			do {
				
				t *= cgp.lsb;
				f1 = Objective(ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data, t, rmse, cgp);
				if (f1 <= f0 - (cgp.lsa * cabs(g0.dotc(dx))))
					break;
				i++;


			} while (i < cgp.lsiter);
			
			printf ("%i, obj: %f, RMS: %f, L-S: %i\n", k, f1, rmse, i);
			
			if (i == cgp.lsiter) {
				printf ("Reached max line search,.... not so good... might have a bug in operators. exiting... \n"); 
				//return;
			}
			
			if (i > 2)
				t0 *= cgp.lsb;
			else if (i < 1)
				t0 /= cgp.lsb;
			
			x += (dx * t);
			
			/* CG computation */
			
			g1  =  Gradient (x, data, cgp);
			bk  =  creal(g1.dotc(g1) / g0.dotc(g0));
			g0  =  g1;
			dx  = -g1; 
			dx += (dx * bk);
			k++;
			
			float dxn = dx.Norm().real();

			if ((k > 8) || dxn < 1e-30) break;
			
		} while (true);
		
	}
	
	
}
#endif /* __COMPRESSED_SENSING_H__ */


/*

#ifdef HAVE_MAT_H	
			MATFile* mf = matOpen ("preobj.mat", "w");
			if (mf == NULL) {
				printf ("Error creating file %s\n", "");
				return;
			}
			x0.MXDump (mf, "x0", ""); ffdbx.MXDump (mf, "ffdbx", ""); ffdbg.MXDump (mf, "ffdbg", ""); 
			ttdbx.MXDump (mf, "ttdbx", ""); ttdbg.MXDump (mf, "ttdbg", "");	g0.MXDump    (mf, "g0",    "");
			if (matClose(mf) != 0) {
				printf ("Error closing file %s\n", "");
				return;
			}
#endif
			
			
*/
