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

#ifndef __DUMMY_RECON_H__
#define __DUMMY_RECON_H__

#include "ReconStrategy.h"

/**
 * @brief Empty recon for test purposes
 */
class DummyRecon : public ReconStrategy {


public:

	/**
	 * @brief Default constructor
	 */
	DummyRecon  () {
		std::cout << "=========== We're actually called =============" << endl; 
	};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~DummyRecon () {};
	
	/**
	 * @brief Do nothing 
	 */
	virtual RRSModule::error_code
	ProcessData () { 

		std::cout << "=========== We're actually called =============" << endl; 

		return RRSModule::OK;

	};
	
};

#endif /* __DUMMY_RECON_H__ */


