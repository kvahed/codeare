/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
 *
 *  Stolen ;) from sparse MRI 0.2 by Michael Lustig
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
#include "FT.hpp"
#include "DWT.hpp"
#include "TVOP.hpp"
#include "CX.hpp"
#include "Algos.hpp"

#include <pthread.h>
/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {


	/**
	 * @brief CG parameters
	 */
	struct CGParam {
		
		int    fft;
		int    lsiter;
		int    lsto ;
		int    cgiter;
		
		double pnorm;
		double tvw;
		double xfmw;
		double cgconv;
		double l1;
		double lsa;
		double lsb;

		DWT*       dwt;
		FT<cxfl>*  ft;
		TVOP*      tvt;
		
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
		int            m_wf;
		int            m_wm;

	};


	float Obj (Matrix<cxfl>& ffdbx, Matrix<cxfl>& ffdbg, Matrix<cxfl>& data, float t) {
	
		Matrix<cxfl> om;
		float        o = 0.0;
		
		om  = ffdbx;
		if (t > 0.0) 
			om += t * ffdbg;
		om -= data;

		o = creal(om.dotc(om));

		return o;

	}


	float ObjTV (Matrix<cxfl>& ttdbx, Matrix<cxfl>& ttdbg, float t, CGParam& cgp) {
		
		Matrix<cxfl> om;
		float        o = 0.0, p = (float)cgp.pnorm/2.0;
		
		om  = ttdbx;
		if (t > 0.0)
			om += t * ttdbg;
		om *= conj(om);
		om += cgp.l1;
		om ^= p;
		
		for (size_t i = 0; i < om.Size(); i++) o += om[i].real();
		
		return cgp.tvw * o;

	}


	/**
	 *
	 */
	float ObjXFM (Matrix<cxfl>& x, Matrix<cxfl>& g, float t, CGParam& cgp) {
		
		Matrix<cxfl> om;
		float        o = 0.0, p = (float)cgp.pnorm/2.0;

		om  = x;
		if (t > 0)
			om += t * g;
		om *= conj(om);
		om += cgp.l1;
		om ^= p;
		
		for (size_t i = 0; i < om.Size(); i++) o += om[i].real();

		return cgp.xfmw * o;

	} 


	float Objective (Matrix<cxfl>& ffdbx, Matrix<cxfl>& ffdbg, Matrix<cxfl>& ttdbx, 
					 Matrix<cxfl>& ttdbg, Matrix<cxfl>&     x, Matrix<cxfl>&     g, 
					 Matrix<cxfl>&  data, float             t, float&         rmse,
					 CGParam&        cgp) {
		
		float obj = 0.0;
		float nz = (float) nnz (data); 
		
		obj  =              Obj    (ffdbx, ffdbg, data, t);
		rmse = sqrt(obj/nz);
		
		obj += (cgp.tvw)  ? ObjTV  (ttdbx, ttdbg, t, cgp) : 0.0;
		obj += (cgp.xfmw) ? ObjXFM (x,     g,     t, cgp) : 0.0;

		return obj;

	}


	/**
	 * @brief Compute gradient of the data consistency
	 */
	Matrix<cxfl> 
	GradObj (Matrix<cxfl>& x, Matrix<cxfl>& data, CGParam& cgp) {
		
		FT<cxfl>& ft = *(cgp.ft);
		DWT& dwt = *(cgp.dwt);

		Matrix<cxfl> g;
		
		g  = ft * (dwt->*x);
		g -= data;
		g  = dwt * (ft->*g);

		return (2.0 * g);

	}


	/**
	 * @brief Compute gradient of L1-transform operator
	 *
	 *
	 */
	Matrix<cxfl> 
	GradXFM   (Matrix<cxfl>& x, CGParam& cgp) {
		
		Matrix<cxfl> g;

		g  = x * conj(x);
		g += cxfl(cgp.l1);
		g ^= (((float)cgp.pnorm)/2.0-1.0);
		g *= x;

		return (cgp.xfmw * g);

	}


	/**
	 * @brief Compute gradient of the total variation operator
	 */
	Matrix<cxfl> 
	GradTV    (Matrix<cxfl>& x, CGParam& cgp) {

		DWT&  dwt = *cgp.dwt;
		TVOP& tvt = *cgp.tvt;
		float p   = ((float)cgp.pnorm)/2.0-1.0;

		Matrix<cxfl> dx, g;

		dx = tvt * (dwt->*x);

		g  = dx * conj(dx);
		g += cxfl(cgp.l1);
		g ^= p;
		g *= dx;
		g *= cxfl(cgp.pnorm);
		g  = dwt * (tvt->*g);

		return (cgp.tvw * g);

	}


	Matrix<cxfl> 
	Gradient (Matrix<cxfl>& x, Matrix<cxfl>& data, CGParam& cgp) {

		Matrix<cxfl> g;
		
		g = GradObj (x, data, cgp);


		if (cgp.xfmw)
			g += GradXFM (x, cgp);

		if (cgp.tvw)
			g += GradTV  (x, cgp);

		return g;

	} 


	void 
	NLCG (Matrix<cxfl>& x, Matrix<cxfl>& data, CGParam& cgp) {


		float        t0 = 1.0, t = 1.0, z = 0.0;
		float        xn = creal(norm(x));
		float        rmse, bk, f0, f1;

		Matrix<cxfl> g0, g1, dx, ffdbx, ffdbg, ttdbx, ttdbg, wx, wdx;

		DWT&      dwt = *cgp.dwt;
		FT<cxfl>& ft  = *cgp.ft;
		TVOP&     tvt = *cgp.tvt;

		g0 = Gradient (x, data, cgp);
		dx = -g0;
		
		for (size_t k = 0; k < cgp.cgiter; k++) {
			
			t = t0;
			
			wx  = dwt->*x;
			wdx = dwt->*dx;

			ffdbx = ft * wx;
			ffdbg = ft * wdx;
			
			if (cgp.tvw) {
				ttdbx = tvt * wx;
				ttdbg = tvt * wdx;
			}
			
			f0 = Objective (ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data, z, rmse, cgp);
		   
			int i = 0;

			while (i < cgp.lsiter) {
				
				t *= cgp.lsb;
				f1 = Objective(ffdbx, ffdbg, ttdbx, ttdbg, x, dx, data, t, rmse, cgp);
				if (f1 <= f0 - (cgp.lsa * t * cabs(g0.dotc(dx))))
					break;
				i++;
				
			} 
			
			printf ("    %02zu - nrms: %1.7f, l-search: %i, ", k, rmse, i); fflush (stdout);

			if (i == cgp.lsiter) {
				printf ("Reached max line search, exiting... \n"); 
				return;
			}
			
			if      (i > 2) t0 *= cgp.lsb;
			else if (i < 1) t0 /= cgp.lsb;

			// Update image
			x += (dx * t);
			
			// CG computation 
			g1  =  Gradient (x, data, cgp);
			bk  =  creal(g1.dotc(g1) / g0.dotc(g0));
			g0  =  g1;
			dx  = -g1 + dx * bk;

			float dxn = creal(norm(dx))/xn;
			printf ("dxnrm: %0.4f\n", dxn);
			if (dxn < cgp.cgconv) 
				break;
			
		} 

	}
	
}
#endif /* __COMPRESSED_SENSING_H__ */

