/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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

#include "nfft3util.h"
#include "nfft3.h"
#include "NCSENSE.hpp"

typedef std::complex<float> raw;

static const int NTHREADS = 8;

namespace RRStrategy {
	
	/**
	 * @brief Conjugate gradient Non-Cartesian SENSE 
	 */
	class CGSENSE : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		CGSENSE () {};
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~CGSENSE ();
		
		/**
		 * @brief Process conjugate gradient SENSE
		 */
		virtual RRSModule::error_code
		Process ();

		/**
		 * @brief Prepare conjugate gradient SENSE
		 */
		virtual RRSModule::error_code
		Prepare ();

		/**
		 * @brief Initialise NuFFT plans
		 */
		virtual RRSModule::error_code
		Init ();
		
		/**
		 * @brief Clean up
		 */
		virtual RRSModule::error_code
		Finalise ();
		
	private:

		NCSENSE<cxfl>* m_ncs;
		
		int                  m_iter;            /**< Maximum number of NuFFT solver iterations */
		int                  m_verbose;         /**< Verbose should give back the reconstruction series? */
		double               m_noise;           /**< Add noise? */
		int                  m_testcase;        /**< Test case. Generate forward data first. */
		double               m_lambda;          /**< Tikhonov factor*/

		double               m_fteps;         /**< NuFFT convergence criterium     */
		int                  m_ftmaxit;           /**< Maximum number of NuFFT solver iterations */
		double               m_cgeps;           /**< CG SENSE convergence criterium  */
		int                  m_cgmaxit;         /**< Maximum number of CG iterations */
		
	};


}
#endif /* __CGSENSE_H__ */
