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

#ifndef __SENSE_HPP__
#define __SENSE_HPP__

#include "ReconStrategy.hpp"
#include "CSENSE.hpp"

/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

	/**
	 * @brief SENSitivity Encoding PPI
	 */
	class SENSE : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		SENSE  () {};
		
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~SENSE () {};
		
		
		/**
		 * @brief Initialise reco
		 */
		virtual RRSModule::error_code
		Init ();
		

		/**
		 * @brief (Re-)prepare module
		 */
		virtual RRSModule::error_code
		Prepare ();

		
		/**
		 * @brief Process data
		 */
		virtual RRSModule::error_code
		Process ();

		
		/**
		 * @brief Clean up
		 */
		virtual RRSModule::error_code
		Finalise ();
		

	private:

		CSENSE<float>*  m_cs;      /**< Cartesian sense operators (Multi-Core Reco)*/

		unsigned short m_ncpus;   /**< Number of threads */
		unsigned short m_af;      /**< Acceleration factor */

		bool           m_compgfm; /**< Compute g-factor map */

	};

}
#endif /* __SENSE_H__ */

