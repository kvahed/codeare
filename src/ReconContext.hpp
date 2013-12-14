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

#ifndef __REMOTE_CONTEXT_HPP__
#define __REMOTE_CONTEXT_HPP__

#include "ReconStrategy.hpp"
#include <vector>


namespace RRStrategy {


	/**
	 * @brief Context of a reconstruction method. Abstraction layer to algorithm backends. 
	 */
	class ReconContext {
		
		
		
	public:
		
		
		/**
		 * @brief        Unload library and destruct.
		 */ 
		~ReconContext    ();
		
		
		/**
		 * @brief        Construct with a strategy name.
		 *               Loads and initialises algorithm. 
		 *               Needs a present config [@see ReconServant::config(const char*)]).
		 *
		 * @name         Name of algorithm
		 */
		ReconContext     (const char* name);
		
		
		/**
		 * @brief        Direct access pointer to underlying algorithm.
		 *
		 * @return       Algorithm
		 */
		inline ReconStrategy*
		Strategy         ();
		
		
		/**
		 * @brief        Process. @see ReconStrategy::Process()
		 *
		 * @return       Success
		 */
		codeare::error_code
		Process          ();
		
		
		/**
		 * @brief        Initialise. @see ReconStrategy::Init()
		 *
		 * @return       Success
		 */
		codeare::error_code
		Init             ();
		
		
		/**
		 * @brief        Prepare. @see ReconStrategy::Prepare()
		 *
		 * @return       Success
		 */
		codeare::error_code
		Prepare          ();
		
		
		/**
		 * @brief        Finalise. @see ReconStrategy::Finalise()
		 *
		 * @return       Success
		 */
		codeare::error_code
		Finalise     ();
		
		
		/**
		 * @brief        @see Reconstrategy::SetConfig(const char*)
		 *
		 * @param  cstr  Serialised XML configuration
		 */
		void
		SetConfig        (const char* cstr);
		
		
		/**
		 * @brief        @see ReconStrategy::ReadConfig(const char*)
		 *
		 * @param  fname File name
		 */
		void
		ReadConfig       (const char* fname);
		
		
		/**
		 * @brief        @see ReconStrategy::SetCXFL(const std::string name, const cxfl_data&)
		 *
		 * @param  name  Name
		 * @param  t     Complex data sequence
		 */
		template <class T> void
		SetMatrix        (const std::string name, const T& t) {

			Workspace::Instance().SetMatrix(name, t);

		}
		
		
		/**
		 * @brief        @see ReconStrategy::GetCXFL(const std::string name, cxfl_data&)
		 *
		 * @param  name  Name
		 * @param  t     Complex data sequence
		 */
		template <class T> void
		GetMatrix        (const std::string name, T& t);
		
		
 		/**
		 * @brief       @see ReconStrategy::Name(const char*)
		 *
		 * @param name  Name
		 */
		void Name (const char* name);
		
		
 		/**
		 * @brief       @see ReconStrategy::Name()
		 *
		 * @return      Name
		 */
		const char* Name ();
		

	private:
		
		
		ReconStrategy*  m_strategy;   /**< Active strategy           */
		void*           m_dlib;       /**< Handle on startegy module */
		
		
		/**
		 * @brief        Default Constructor.
		 */
		ReconContext     ();
		
		
	};
	
}
#endif 
