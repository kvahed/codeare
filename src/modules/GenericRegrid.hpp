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

#ifndef __GENERIC_REGRID_HPP__
#define __GENERIC_REGRID_HPP__

#include "ReconStrategy.hpp"

namespace RRStrategy {
	
	/**
	 * @brief Generic regridder
	 *        Expects to identically sized matrices (data & k-space positions)
	 */
	class GenericRegrid : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		GenericRegrid  () {};
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~GenericRegrid () {};
		
		/**
		 * @brief Regrid data to Cartesian k-space
		 */
		virtual error_code
		Process () {
			return OK;
		};
		
		/**
		 * @brief Do nothing 
		 */
		virtual error_code
		Init () {
			return OK;
		}

		/**
		 * @brief Do nothing 
		 */
		virtual error_code
		Finalise () {
			return OK;
		}
		
	};
	
}

#endif /* __GENERIC_REGRID_HPP__ */
