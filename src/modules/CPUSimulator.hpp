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
		~CPUSimulator ();
		
		
		/**
		 * @brief       Default constructor
		 */
		CPUSimulator  (SimulationBundle* sb);


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
		Simulate  ();

		virtual void 
		Simulate  (const bool& mode);
		

	protected:
		
		Matrix<cplx> m_sig;
		double       m_gdt;
		size_t       m_nt;
		size_t       m_nc;
		
	};
		
}

#endif // __CPU_SIMULATOR_HPP__
