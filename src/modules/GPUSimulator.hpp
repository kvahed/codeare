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

#ifndef __GPU_SIMULATOR_HPP__
#define __GPU_SIMULATOR_HPP__

#include "SimulationStrategy.hpp"

namespace RRStrategy {
	
	/**
	 * @brief Base class for simulation stratgies used by direct method
 	 */
	class GPUSimulator : public SimulationStrategy {
		
		
	public:
		
		
		/**
		 * @brief       Default destructor
		 */
		virtual 
		~GPUSimulator () {};
		
		
		/**
		 * @brief       Default constructor
		 */
		GPUSimulator  () {};


		/**
		 * @brief       Piece-wise constant bloch simulation
		 *
		 * @see         SimulationStrategy::Simulate()
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
		 * @see         SimulationStrategy::Simulate()
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
		 * @see         SimulationStrategy::Simulate()
		 */
		virtual void
		SimulateExc    (const Matrix<cplx>&   txm, const Matrix<cplx>&   rf, 
						const Matrix<double>&  gr, const Matrix<double>&  r, 
						const Matrix<double>& b0m, const double&         dt, 
						const bool&            v,  const size_t&        pos, 
						const int&            tid,       Matrix<double>&  m);
		
		
	};
		
}

#endif // __GPU_SIMULATOR_HPP__
