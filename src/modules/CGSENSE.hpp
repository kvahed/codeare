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

#ifndef __CGSENSE_HPP__
#define __CGSENSE_HPP__

#include "ReconStrategy.hpp"
#include "Algos.hpp"

#include "NCSENSE.hpp"

typedef std::complex<float> raw;

static const int NTHREADS = 4;

namespace RRStrategy {
	
	/**
	 * @brief Conjugate gradient Non-Cartesian SENSE 
	 */
	class CGSENSE : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		CGSENSE () : m_cgeps(1.0e-7), m_fteps(1.0e-3), m_cgmaxit(10),
					 m_ftmaxit(3), m_noise(0.0), m_lambda(1.0e-6), m_testcase(0),
					 m_verbose(0), m_nthreads (0), m_m(1), m_nk(1) {}
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~CGSENSE ();
		
		/**
		 * @brief Process conjugate gradient SENSE
		 */
		virtual codeare::error_code
		Process ();

		/**
		 * @brief Prepare conjugate gradient SENSE
		 */
		virtual codeare::error_code
		Prepare ();

		/**
		 * @brief Initialise NuFFT plans
		 */
		virtual codeare::error_code
		Init ();
		
		/**
		 * @brief Clean up
		 */
		virtual codeare::error_code
		Finalise ();
		
	private:

		NCSENSE<cxfl>  m_ncs;
		
		int             m_verbose;   /**< Verbose should give back the reconstruction series? */
		int             m_testcase;  /**< Test case. Generate forward data first.             */
		int             m_ftmaxit;   /**< Maximum number of NuFFT solver iterations           */
		int             m_cgmaxit;   /**< Maximum number of CG iterations                     */
		int             m_nthreads;  /**< Number of threads                                   */
		int             m_nk;        /**< Number of kspace samples                            */
        int             m_m;

        bool            m_3rd_dim_cart; /**< 3rd NUFFT direction is Cartesian (stack(spirals/stars)) */
		
		double          m_noise;     /**< Add noise?                                          */
		double          m_lambda;    /**< Tikhonov factor                                     */
		double          m_fteps;     /**< NuFFT convergence criterium                         */
		double          m_cgeps;     /**< CG SENSE convergence criterium                      */

	};


}
#endif /* __CGSENSE_H__ */
