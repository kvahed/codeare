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

#define __CL_ENABLE_EXCEPTIONS

#include "cl.hpp"

#define MAXDEVS 1

namespace RRStrategy {
	

	/**
	 * @brief Simple time equidistant bloch simulator for OpenCL 
 	 */
	class GPUSimulator : public SimulationStrategy {
		
		
	public:
		
		
		/**
		 * @brief       Default destructor
		 */
		virtual 
		~GPUSimulator  () {};
		
		
		/**
		 * @brief       Default constructor
		 */
		GPUSimulator   (SimulationBundle* sb);


		/**
		 * @brief        Piece-wise constant bloch simulation
		 *
		 * @see          SimulationStrategy::Simulate()
		 */
		virtual void 
		Simulate        ();
		
		
	protected:
		

		/**
		 * @brief        Prepare GPU processing<br/>(i.e. Load kernel & transfer data to GPU memory)
		 */
		bool
		PrepareKernel   ();		
		

		/**
		 * @brief        Convenience for human readable error
		 *
		 * @param   e    Error code
		 * @return       Human readable error string
		 */
		const char*	
		ErrorString     (cl_int e);		


		/**
		 * @brief        Read OpenCL program source from file
		 *
		 * @param  fname File name
		 * @param  size  File size
		 * @return       File content
		 */
		char* 
		ReadSource      (const char* fname, int* size);		


		/**
		 *  @brief       Build OpenCL program
		 *
		 *  @param  ksrc Kernel source
		 */
		bool
		BuildProgram    (std::string ksrc);
		

		/**
		 *  @brief       Allocate memory on device and copy data
		 */
		void 
		SetDeviceData   ();
		

		/**
		 *  @brief       Run kernel
		 */
		void 
		RunKernel       ();
		

		/**
		 *  @brief       Get data back from device and clean up
		 */
		void 
		GetDeviceData   ();
		

        unsigned int              m_dev;    /**!<   */
		std::vector<cl::Device>   m_devs;   /**!<   */
        cl::Context               m_ctxt;   /**!<   */
        cl::CommandQueue          m_cmdq;   /**!<   */
        cl::Program               m_prg;    /**!<   */
        cl_int                    m_error;  /**!<   */
        cl::Event                 m_event;  /**!<   */
        cl::Kernel                m_kernel; /**!<   */

		// RO device
		cl::Buffer                ocl_tb1;
		cl::Buffer                ocl_sb1;
		cl::Buffer                ocl_agr;
		cl::Buffer                ocl_tr;
		cl::Buffer                ocl_sr;
		cl::Buffer                ocl_tb0;
		cl::Buffer                ocl_sb0;
		cl::Buffer                ocl_tm;
		cl::Buffer                ocl_sm;
		cl::Buffer                ocl_jac;
		cl::Buffer                ocl_gdt;
		cl::Buffer                ocl_n;

		// RW on device
		cl::Buffer                ocl_rf;
		cl::Buffer                ocl_m;

		float                     m_gdt;
		size_t                    m_nt;
		size_t                    m_nc;
		size_t                    m_na;
		size_t                    m_ne;
		size_t                    m_nl;


	};

	
	
	
}

#endif // __GPU_SIMULATOR_HPP__
