/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Daniel Brenner
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

#ifndef __AFI_HPP__
#define __AFI_HPP__

#include "ReconStrategy.hpp"

using namespace RRServer;


/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

	/**
	 * @brief AFI reconstruction
	 */
	class AFI : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		AFI  () : 
			m_use_real (true),
			m_retain_phase (true) {};
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~AFI () {};
		
		
		/**
		 * @brief Do nothing 
		 */
		virtual RRSModule::error_code
		Process ();

		
		/**
		 * @brief Do nothing 
		 */
		virtual RRSModule::error_code
		Init ();
		

		/**
		 * @brief Do nothing 
		 */
		virtual RRSModule::error_code
		Finalise ();


	private:

		bool m_use_real;
		bool m_retain_phase;
		
	};

}
#endif /* __AFI_H__ */

