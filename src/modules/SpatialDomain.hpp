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

#ifndef __SPATIAL_DOMAIN_HPP__
#define __SPATIAL_DOMAIN_HPP__

#include "ReconStrategy.hpp"

using namespace RRServer;


/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

	/**
	 * @brief Empty recon for test purposes
	 */
	class SpatialDomain : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		SpatialDomain  ();
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~SpatialDomain ();
		
		
		/**
		 * @brief Process
		 */
		virtual RRSModule::error_code
		Process ();

		
		/**
		 * @brief Initialise
		 */
		virtual RRSModule::error_code
		Init ();
		

		/**
		 * @brief Finalise
		 */
		virtual RRSModule::error_code
		Finalise ();



	private:

		int    m_nc;     /**<Transmit channels    */
		int*   m_pd;     /**< Pulse durations     */
		int    m_gd;     /**< Gradient durations  */
		int    m_ns;     /**< # Spatial positions */
		int    m_nk;     /**< # kt-points         */
		int    m_maxiter; /**< # Variable exchange method iterations */

		double m_lambda; /**< Tikhonov parameter  */
		double m_rflim;

	};

}
#endif /* __DUMMY_RECON_H__ */


