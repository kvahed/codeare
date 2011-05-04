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

#ifndef __RECONCONTEXT_HPP__
#define __RECONCONTEXT_HPP__

#include <vector>

#include "ReconStrategy.hpp"

namespace RRServer {

	/**
	 * @brief Context of a reconstruction method
	 */
	class ReconContext {
		
		
		
	public:
		
		
		/**
		 * @brief Default Constructor
		 */
		ReconContext () {}
		
		
		/**
		 * @brief Invoce destruction on my startegy and exit
		 */ 
		~ReconContext ();
		
		
		/**
		 * @brief Construct with a strategy
		 */
		ReconContext (const char* name);
		
		
		/**
		 * @brief get active startegy
		 */
		inline ReconStrategy*
			Strategy     () {
			return m_strategy;
		}
		
		/**
		 * @brief Process data with given strategy
		 */
		RRSModule::error_code
			Process () {
			return m_strategy->Process();
		}
		
		
	private:
		
		ReconStrategy*            m_strategy;   /**< Active strategy           */
		void*                     m_dlib;       /**< Handle on startegy module */
		
	};
	
}
#endif 
