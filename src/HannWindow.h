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

#ifndef __HANN_WINDOW_H__
#define __HANN_WINDOW_H__

#include "ReconStrategy.h"

/**
 * @brief Hann window
 */
class HannWindow : public ReconStrategy {


public:

	/**
	 * @brief Default constructor
	 */
	HannWindow  () {};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~HannWindow () {};
	
	/**
	 * @brief Apply Hann window to k-space
	 */
	virtual RRSModule::error_code
	ProcessData () {
		return RRSModule::OK;
	};
	
};

#endif /* __HANN_WINDOW_H__ */

extern "C" {

	ReconStrategy* maker() {
		return new HannWindow;
	}

	class Proxy { 

	public:
		
		Proxy(){
			// register the maker with the factory
			factory["HannWindow"] = maker;
		}
	};
	
	// our one instance of the proxy
	Proxy p;
}
