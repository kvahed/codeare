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

#ifndef __RECONCONTEXT_H__
#define __RECONCONTEXT_H__

#include <vector>
#include <dlfcn.h>
#include "ReconStrategy.h"

using namespace std;

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
	~ReconContext () {
		
		// load destroy symbol
		destroy_t* destroy = (destroy_t*) dlsym(m_dlib, "destroy");
		err                = dlerror();

		if (err)
			cerr << "Cannot load symbol destroy: " << err << '\n';
		
		destroy(m_strategy);
		dlclose(m_dlib);

	};


	/**
	 * @brief Construct with a strategy
	 */
	ReconContext (const char* name) {

		m_dlib = dlopen(name, RTLD_NOW);

		err = dlerror();

		if (!m_dlib) 
			cerr << "Cannot load library: " << err << '\n';
		
		// reset errors
		err = dlerror();
		
		// load create symbol
		create_t* create = (create_t*) dlsym(m_dlib, "create");

		err = dlerror();

		if (err) 
			cerr << "Cannot load symbol create: " << err << '\n';
		
		m_strategy = create();

	}
	
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
	inline RRSModule::error_code
	ProcessData () {
		return m_strategy->ProcessData();
	}
	

private:

	ReconStrategy*            m_strategy;   /**< Active strategy      */
	void*                     m_dlib;
	
	const char*               err; 

};

#endif 
