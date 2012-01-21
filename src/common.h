/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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

#ifndef __COMMON_H__
#define __COMMON_H__

/**
 * Returned error codes 
 */
enum error_code   {
	
	OK,
	ERROR_GENERAL,
	ERROR_UNKNOWN_METHOD,
	FATAL_SYSTEM_CALL,
	LAPACK_WORKSPACE_QUERY_FAILED,
	CG_DID_NOT_CONVERGE,
	CANNOT_LOAD_LIBRARY,
	UNSUPPORTED_DIMENSION,
	UNSUPPORTED_IMAGE_MATRIX,
	ZERO_NODES,
	CGSENSE_ZERO_CHANNELS,
	FILE_ACCESS_FAILED,
	CONTEXT_CONFIGURATION_FAILED,
	CONTEXT_NOT_FOUND

};

enum coords {
	X,
	Y,
	Z
};

enum sim_mode {
	ACQUIRE,
	EXCITE
};

#endif //__COMMON_H__
