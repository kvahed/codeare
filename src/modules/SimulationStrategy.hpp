/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
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

#ifndef __SIMULATION_STRATEGY_HPP__
#define __SIMULATION_STRATEGY_HPP__

enum SIM_MODE {
	AO,  // Acquisition only
	AE,  // Acquire and excite
	CG   // CG ATA
};

#include "common.h"
#include "Matrix.hpp"
#include "IOContext.hpp"

namespace RRStrategy {
	
	/**
	 * @brief  Simulation package structure
	 */
	struct SimulationBundle {

		// Incoming
		boost::shared_ptr<Matrix<cxfl> >    b1;  /**<! b1                               */ 
		
		boost::shared_ptr<Matrix<float> >    g;  /**<! Acquisition gradients            */ 

		boost::shared_ptr<Matrix<float> >    r;  /**<! spatial vectors                  */ 

		boost::shared_ptr<Matrix<float> >   b0;  /**<! b0 maps                          */ 

		boost::shared_ptr<Matrix<cxfl> >  tmxy;  /**<! starting magnetisation (target)  */
		boost::shared_ptr<Matrix<float> >  tmz;
		boost::shared_ptr<Matrix<cxfl> >  smxy;  /**<!                        (sample)  */
		boost::shared_ptr<Matrix<float> >  smz;  
		boost::shared_ptr<Matrix<float> >  roi;  /**<! ROI                    (sample)  */

		boost::shared_ptr<Matrix<float> >  jac;  /**<! jacobian j(k(t))                 */

		int                    np;  /**<! # threads                        */
		int                  mode;  /**<! mode (0:single run, 1:iterative) */

		float                  dt;  /**<! time step                        */
		float               cgeps;  /**<! CGNR convergence criterium       */
		float              lambda;  /**<! Tikhonov regularisation factor   */
		int                  cgit;  /**<! CGNR iterations                  */
		
		bool                    v;  /**<! verbose                          */
		bool                  cb0;  /**<! correct b0?                      */
		
		// Outgoing
		boost::shared_ptr<Matrix<cxfl> >    rf;  /**<! RF pulses                         */
		boost::shared_ptr<Matrix<cxfl> >   mxy;  /**<! Excited transverse magnetisation  */
		boost::shared_ptr<Matrix<float> >   mz;  /**<! Longitudinal magnetisation        */

	};

	/**
	 * @brief Base class for simulation stratgies used by direct method
 	 */
	class SimulationStrategy {		
		
	public:
		

		/**
		 * @brief       Construct with bundle
		 */
		SimulationStrategy  (SimulationBundle* sb) { m_sb = sb; }

		/**
		 * @brief       Default destructor
		 */
		virtual 
		~SimulationStrategy ()                     {};


		/**
		 * @brief       Piece-wise constant bloch simulation
		 */
		virtual void 
		Simulate  () = 0;

		
	protected: 
		
		SimulationStrategy  () : m_sb(0)                   {};

		SimulationBundle* m_sb;

	};
		




}

#endif // SimulationStrategy
