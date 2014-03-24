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

#ifndef __COMMON_H__
#define __COMMON_H__

#include <cstddef>

/**
 * Returned error codes 
 */
namespace codeare {
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
	CONTEXT_NOT_FOUND,
	MEM_ALLOC_FAILED,
	FILE_OPEN_FAILED,
	FILE_READ_FAILED,
	FILE_WRITE_FAILED,
	NULL_FILE_HANDLE,
	UNIMPLEMENTED_METHOD,
	HDF5_ERROR_FOPEN,
	HDF5_ERROR_FCLOSE,
	HDF5_ERROR_FFLUSH,
	GENERAL_IO_ERROR,
    CONFIG_LOAD_FILE_FAILED,
    WRONG_OR_NO_DATASET,
    CONFIG_MISSING_CHAIN,
    CONFIG_EMPTY_CHAIN,
    NULL_STRATEGY,
    NO_MATRIX_IN_WORKSPACE_BY_NAME,
    WRONG_MATRIX_TYPE

};
}
enum coords {
	X, Y, Z
};

enum sim_mode {
	ACQUIRE, EXCITE
};

enum connection_type {
	CON_LOCAL, CON_REMOTE
};

enum data_type {
	T_CXFL,
	T_CXDB,
	T_FLOAT,
	T_DOUBLE,
	T_SHORT,
	T_LONG,
	T_BOOL
};



/**
 * Some constants
 */
// PI
#ifndef PI
    # define PI 3.141592653589793238462643383279502884197169399375105820974944592
#endif

#ifndef TWOPI
    #define TWOPI 6.283185307179586476925286766559005768394338798750211641949889185
#endif

// Gamma in Hz
#ifndef GAMMA
    #define GAMMA 4.257748e1
#endif

// Gamma in radians
#ifndef RGAMMA
    #define RGAMMA 2.67513e2
#endif

// Bytes & Co
#define KB          1024;
#define MB       1048576;
#define GB    1073741824;
#define TB 1099511627776;

#define IZERO 0
#define IONE  1

#if defined (_MSC_VER) && _MSC_VER<1300
template<class T> inline static T
c_multiply (const T a, const T b);

template<> inline static size_t
c_multiply (const size_t a, const size_t b) {
    return (a*b);
}
#else
template<class T> inline static T
c_multiply (const T a, const T b) {
    return (a*b);
}
#endif

#if defined(_MSC_VER)
  #define JL_SIZE_T_SPECIFIER    "%Iu"
  #define JL_SSIZE_T_SPECIFIER   "%Id"
  #define JL_PTRDIFF_T_SPECIFIER "%Id"
#elif defined (__GNUC__)
  #define JL_SIZE_T_SPECIFIER    "%zu"
  #define JL_SSIZE_T_SPECIFIER   "%zd"
  #define JL_PTRDIFF_T_SPECIFIER "%zd"
#else
// TODO figure out which to use.
  #if NUMBITS == 32
    #define JL_SIZE_T_SPECIFIER    something_unsigned
    #define JL_SSIZE_T_SPECIFIER   something_signed
    #define JL_PTRDIFF_T_SPECIFIER something_signed
  #else
    #define JL_SIZE_T_SPECIFIER    something_bigger_unsigned
    #define JL_SSIZE_T_SPECIFIER   something_bigger_signed
    #define JL_PTRDIFF_T_SPECIFIER something-bigger_signed
  #endif
#endif

#endif //__COMMON_H__
