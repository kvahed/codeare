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

#ifndef __CPU_SIMULATOR_HPP__
#define __CPU_SIMULATOR_HPP__

#include "SimulationStrategy.hpp"

namespace RRStrategy {
	
	/**
	 * @brief Simple time equidistant bloch simulator 
 	 */
	class CPUSimulator : public SimulationStrategy {
		
		
	public:
		
		
		/**
		 * @brief       Clean up and destroy
		 */
		virtual 
		~CPUSimulator  ();
		
		
		/**
		 * @brief       Construct with bundle
		 *
		 * @param  sb   Sumulation bundle
		 */
		CPUSimulator   (SimulationBundle* sb);


		/**
		 * @brief       Piece-wise constant bloch simulation
		 */
		virtual void 
		Simulate       ();



	protected:
		

		float          m_gdt; /*<! \gamma*dt          */
		float          m_rfsc; /*<! RF scale */
		size_t         m_nt;  /*<! # timepoints       */
		size_t         m_nc;  /*<! # channels         */


	};
		
}

#endif // __CPU_SIMULATOR_HPP__
