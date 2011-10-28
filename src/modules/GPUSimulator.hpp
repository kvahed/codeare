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
#include "../config.h"

#include <vector>

#if defined HAVE_CL_CL_H
#define NVIDIA
#elif defined HAVE_OPENCL_CL_H
#define APPLE
#endif

#include "cl.hpp"

#define __CL_ENABLE_EXCEPTIONS


namespace RRStrategy {
	
	const unsigned int MAX_GPU_COUNT = 8;


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
		GPUSimulator  ();


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
		

		cl_context       m_ctxt;
		cl_kernel        m_simkern[MAX_GPU_COUNT];
		cl_command_queue m_cq[MAX_GPU_COUNT];
		cl_platform_id   m_cpi;
		cl_uint          m_devs;
		cl_device_id*    m_devids;
		cl_int           m_err;

	};
		
}

#endif // __GPU_SIMULATOR_HPP__
