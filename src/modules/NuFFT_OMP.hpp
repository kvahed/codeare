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

#ifndef __NUFFT_OMP_HPP__
#define __NUFFT_OMP_HPP__

#include "nfftstub.h"

#include "ReconStrategy.hpp"

static const int NTHREADS = 8;

using namespace RRServer;

namespace RRStrategy {
	
	/**
	 * @brief Non uniform FFT
	 */
	class NuFFT_OMP : public ReconStrategy {
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		NuFFT_OMP  ();
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~NuFFT_OMP ();
		
		/**
		 * @brief Dump data to disk
		 */
		virtual RRSModule::error_code
		Process ();
		
		/**
		 * @brief Dump data to disk
		 */
		virtual RRSModule::error_code
		Init ();
		
	private:
		
		// We only support 2 dims for the time being
		int m_dim;
		
		// Some more stuff
		int*      m_N;                          /**< @brief Image matrix side length */
		int*      m_n;                          /**< @brief Oversampling */
		int       m_M;                          /**< @brief Number of k-space knots */
		int       m_maxit;                      /**< @brief Number of Recon iterations (NFFT 3) */
		int       m_shots;                      /**< @brief Number of shots */
		int       m_verbose;                    /**< @brief Store all shots separately? */

		double    m_epsilon;                    /**< @brief Convergence criterium */

		nfft_plan           m_fplan[NTHREADS];  /**< NuFFT plan                      */
		solver_plan_complex m_iplan[NTHREADS];  /**< iNuFFT plan                     */
		
	};
	
}
#endif /* __NUFFT_HPP__ */
