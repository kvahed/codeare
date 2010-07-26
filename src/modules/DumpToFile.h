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

#ifndef __DUMP_TO_FILE_H__
#define __DUMP_TO_FILE_H__

#include "../ReconStrategy.h"

/**
 * @brief Dump data to file
 */
class DumpToFile : public ReconStrategy {


public:
	
	/**
	 * @brief Default constructor
	 */
	DumpToFile  () {};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~DumpToFile () {};
	
	/**
	 * @brief Dump data to disk
	 */
	virtual RRSModule::error_code
	ProcessData ();
	
};

#endif /* __DUMMY_RECON_H__ */
