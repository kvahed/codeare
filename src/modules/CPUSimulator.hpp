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
	 * @brief Base class for simulation stratgies used by direct method
 	 */
	class CPUSimulator : public SimulationStrategy {
		
		
	public:
		
		
		/**
		 * @brief       Default destructor
		 */
		virtual 
		~CPUSimulator () {};
		
		
		/**
		 * @brief       Default constructor
		 */
		CPUSimulator  () {};


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
		Simulate  (const Matrix<cplx>&   txm, const Matrix<cplx>&   rxm, 
				   const Matrix<cplx>&    rf, const Matrix<double>&  gr, 
				   const Matrix<double>&   r, const Matrix<double>&  m0, 
				   const Matrix<double>& b0m, const double&          dt, 
				   const bool&           exc, const bool&             v, 
				   const size_t&          np, 
				         Matrix<cplx>&   res, Matrix<double>&         m);
		
		
	private:
		

		/**
		 * @brief       Simulate single shot reception of freely precessing isochromat along gradient trajectory<br/>
		 *              (i.e. forward Fourier transform incl. effect of Receive and b0 maps)<br/>
		 *              Expects res to carry the correct size and dimensions
		 *
		 * INPUT:
		 * @param  rxm  Receive sensitivities
		 * @param  gr   Gradient trajectory
		 * @param  r    Positions
		 * @param  m0   Initial magnetisation state
		 * @param  b0m  B0 map
		 * @param  dt   Timestep
		 * @param  v    Verbose
		 * @param  pos  Position of isochromat in r
		 * @param  tid  Thread ID
		 *
		 * OUTPUT:
		 * @param  res  Resulting signal (output) 
		 */
		virtual void
		SimulateRecv   (const Matrix<cplx>&   rxm, const Matrix<double>& gr, 
						const Matrix<double>&   r, const Matrix<double>& m0, 
						const Matrix<double>& b0m, const double&         dt, 
						const bool&             v, const size_t&        pos, 
						const int&            tid,       Matrix<cplx>&  res);
		
		/**
		 * @brief       Simulate single shot excitation of a single isochromat of r matrix along gradient trajectory<br/>
		 *              (i.e. inverse Fourier transform incl. effect of Receive and b0 maps)<br/>
		 *              Expects res to carry the correct size and dimensions
		 *
		 * INPUT:
		 * @param  txm  Receive sensitivities
		 * @param  rf   Complex RF field
		 * @param  gr   Gradient trajectory
		 * @param  r    Positions
		 * @param  b0m  B0 map
		 * @param  dt   Timestep
		 * @param  v    Verbose
		 * @param  pos  Position of isochromat in r
		 * @param  tid  Thread id
		 *
		 * OUTPUT:
		 * @param  m    Resulting magntisation (output)
		 */
		virtual void
		SimulateExc    (const Matrix<cplx>&   txm, const Matrix<cplx>&   rf, 
						const Matrix<double>&  gr, const Matrix<double>&  r, 
						const Matrix<double>& b0m, const double&         dt, 
						const bool&            v,  const size_t&        pos, 
						const int&            tid,       Matrix<double>&  m);
		
		
	};
		
}

#endif // __CPU_SIMULATOR_HPP__
