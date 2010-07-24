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

#ifndef __INVERT_ORDER_H__
#define __INVERT_ORDER_H__

#include "ReconStrategy.h"

/**
 * @brief Reorder for test purposes
 */
class InvertOrder : public ReconStrategy {


public:

	/**
	 * @brief Default constructor
	 */
	InvertOrder  () {};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~InvertOrder () {};
	
	/**
	 * @brief Invert the order of the data. 
	 *        No good for anything besides testing.
	 */
	virtual RRSModule::error_code
	ProcessData ();
	
};

#endif /* __INVERT_ORDER_H__ */

extern "C" {

	ReconStrategy* maker() {
		return new InvertOrder;
	}

	class Proxy { 

	public:
		
		Proxy(){
			// register the maker with the factory
			//factory["InvertOrder"] = maker;
			m_maker = maker;
		}
	};
	
	// our one instance of the proxy
	Proxy p;
}
