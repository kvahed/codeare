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

#ifndef __FUNCTOR_CONTAINER_H__
#define __FUNCTOR_CONTAINER_H__

#include "ReconContext.hpp"

using namespace RRStrategy;


/**
 * @brief Container for data and reconstructions
 */
class FunctorContainer {



public:


	
	/**
	 * @brief      Default constructor
	 */
	FunctorContainer () {}


	/**
	 * @brief      Destructor
	 */
	virtual
	~FunctorContainer ();

	
	/**
	 * @brief      Process startegy (Needs initialisation @see Init)
	 *
	 * @param name Name of library
	 * @return     Sucess
	 */
	virtual error_code
	Process        (const char* name);
	
	
	/**
	 * @brief      Prepare startegy (Needs initialisation @see Init)
	 *
	 * @param name Name of library
	 * @return     Sucess
	 */
	virtual error_code
	Prepare        (const char* name);
	
	
	/**
	 * @brief      Initialise strategy (Configuration document needs to be set first @see config)
	 * 
	 * @param name Name of processing library
	 * @return     success
	 */
	virtual error_code
	Init          (const char* name);
	
	
	/**
	 * @brief      Finalise algorithm
	 *
	 * @param name Name of processing library
	 */
	virtual error_code
	Finalise       (const char* name = 0);
	
	
	/**
	 * @brief      Clean up left over objects
	 *
	 * @return     Success
	 */
	virtual error_code
	CleanUp        ();
	

	/**
	 * @brief      Set configuration
	 *
	 * @param c    Configuration string
	 */
	virtual void 
	config         (const char* c);
	

	
protected:


	char*                                m_config;   /**< Serialised XML document  */
	std::map<std::string, ReconContext*> m_contexts; /**< Reconstruction contexts (Abstraction layer to algorithms)*/

};

#endif //__FUNCTOR_CONTAINER_H__
