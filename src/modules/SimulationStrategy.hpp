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

#ifndef __SIMULATION_STRATEGY_HPP__
#define __SIMULATION_STRATEGY_HPP__

#include "../common.h"
#include "Matrix.hpp"

namespace RRStrategy {
	
	struct SimulationBundle {

		// Incoming
		Ptr< Matrix<cplx> >   tb1;  /**<! b1                   (target)   */ 
		Ptr< Matrix<cplx> >   sb1;  /**<!                      (sample)   */ 
		
		Ptr< Matrix<double> > agr;  /**<! Acquisition gradients           */ 
		
		Ptr< Matrix<double> > tr;   /**<! spatial vectors        (target) */ 
		Ptr< Matrix<double> > sr;   /**<! */ 

		Ptr< Matrix<double> > tb0;  /**<! b0 maps                (target) */ 
		Ptr< Matrix<double> > sb0;  /**<!                        (sample) */ 

		Ptr< Matrix<double> > tm;   /**<! starting magnetisation (target) */
		Ptr< Matrix<double> > sm;   /**<!                        (sample) */

		Ptr< Matrix<double> > jac;  /**<! jacobian j(k(t))                */

		int                   np;   /**<! # threads                       */
		int                   mode; /**<! mode                            */

		double                dt;   /**<! time step                       */
		
		bool                  v;    /**<! verbose                         */
		
		// Outgoing
		Ptr< Matrix<cplx> >   rf;   /**<! RF pulses                       */
		Ptr< Matrix<double> > magn; /**<! Excited magnetisation           */

		bool Dump (std::string odf) {
			
#ifdef HAVE_MAT_H	
			
			MATFile* mf = matOpen (odf.c_str(), "w");
			
			if (mf == NULL) {
				printf ("Error creating file %s\n", odf.c_str());
				return false;
			}
			
			sb1->MXDump (mf, "sb1");	
			tb1->MXDump (mf, "tb1");
			tr->MXDump (mf, "tr");
			sr->MXDump (mf, "sr");
			tb0->MXDump (mf, "tb0");
			sb0->MXDump (mf, "sb0");
			tm->MXDump (mf, "tm");
			sm->MXDump (mf, "sm");
			jac->MXDump (mf, "jac");
			rf->MXDump (mf, "rf");
			magn->MXDump (mf, "magn");
			
			if (matClose(mf) != 0) {
				printf ("Error closing file %s\n", odf.c_str());
				return false;
			}
#endif
		}
		
	};

	/**
	 * @brief Base class for simulation stratgies used by direct method
 	 */
	class SimulationStrategy {		
		
	public:
		

		/**
		 * @brief       Construct with bundle
		 */
		SimulationStrategy  (SimulationBundle* sb) {m_sb = sb;}

		/**
		 * @brief       Default destructor
		 */
		virtual 
		~SimulationStrategy ()                    { 	printf ("    SimulationStrategy: Cleaned up.\n"); };


		/**
		 * @brief       Piece-wise constant bloch simulation
		 *
		 * INPUT:
		 * @param  txm  Transmit sensitivity    (Nr  x Ntxc)
		 * @param  rxm  Receive sensitivity     (Nr  x Nrxc) 
		 * @param  rf   RF field                (Nt  x Nt  )
		 * @param  gr   Gradient                (1-3 x Nt  )
		 * @param  r    Spatial positions       (1-3 x Nr  )
		 * @param  m0   Starting magnetisation  (3   x Nr  )
		 * @param  b0m  B0 map                  (Nr        )
		 * @param  dt   time step
		 * @param  exc  Exciting or receiving  
		 * @param  v    Verbose                 (Scalar: false = only end, true = all time points)
		 * @param  np   # parallel processes    (scalar)
		 *
		 * OUTPUT:
		 * @param  res  Result of simulation    (Nr  x 3 (x Nt))
		 * @param  m    Resulting magnetisation (3   x Nt)
		 */
		virtual void 
		Simulate  () = 0;

		
	protected: 
		
		SimulationStrategy  ()                    {};

		SimulationBundle* m_sb;

	};
		




}

#endif // SimulationStrategy
