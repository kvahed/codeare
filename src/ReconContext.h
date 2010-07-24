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

#include <dlfcn.h>

#include "ReconStrategy.h"

using namespace std;

/**
 * @brief Context of a reconstruction method
 */
class ReconContext {
  

  
public:


	/**
	 * @brief Construct with setting up available strategies
	 */
	ReconContext ();


	/**
	 * @brief Default destructor
	 */
	~ReconContext () {
		
		delete m_strategy;

		if (m_dlib != NULL)
			dlclose (m_dlib);

	};


	/**
	 * @brief Construct with a strategy
	 */
	ReconContext (ReconStrategy* strategy) {};


	/**
	 * @brief Alter strategy
	 */
	inline void 
	Strategy     (ReconStrategy* strategy) {
		m_strategy = strategy;
	}
	

	/**
	 * @brief Alter strategy
	 */
	inline void 
	Strategy     (method m) {

		//m_strategy = m_strategies.at(m);

	}


	/**
	 * @brief Alter strategy
	 */
	bool
	Strategy     (const char* name) {

		m_dlib = dlopen(name, RTLD_NOW);

		if(m_dlib == NULL){
			cerr << dlerror() << endl;
		}

		void*  maker = dlsym(m_dlib, "maker");
		m_strategy   = static_cast<ReconStrategy *()>(maker)();
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

	ReconStrategy*             m_strategy;    /**< Active strategy      */
	//vector < ReconStrategy* >  m_strategies;  /**< Available strategies */

	void*                      m_dlib;
	map <string, maker_t *, less<string> >           factory;
	map <string, maker_t *, less<string> >::iterator factitr;	
};

#endif 
