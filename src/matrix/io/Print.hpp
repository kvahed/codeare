/*
 *  codeare Copyright (C) 2007-2010 
 *                        Kaveh Vahedipour
 *                        Forschungszentrum Juelich, Germany
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

#ifndef __PRINT_HPP__
#define __PRINT_HPP__

#include "Matrix.hpp"

#include <sstream>
#include <boost/format.hpp>



/**
 * @brief          Print dimensions to string
 * 
 * @param  M       Matrix
 * @return         String of dimensions
 */
template <class T> inline static std::string 
DimsToString (const Matrix<T>& M) {
	std::stringstream ss;
	for (size_t i = 0; i < ndims(M); i++)
		ss << (int)M.Dim(i) << " ";
	return ss.str();
}
	
	
/**
 * @brief          Print dimensions to c string
 *
 * @param  M       Matrix
 * @return         C string of dimensions
 */
template <class T> inline static const char* 
DimsToCString (const Matrix<T>& M) {
	return DimsToString(M).c_str();
}
	
	
/**
 * @brief          Print resolutions to string
 * 
 * @param  M       Matrix
 * @return         String of resolutions
 */
template <class T> inline static std::string 
ResToString (const Matrix<T>& M) {
	std::stringstream ss;
	for (size_t i = 0; i < ndims(M); i++)
		ss << M.Res(i) << " ";
	return ss.str();
}



/**
 * @brief          Print dimensions to c string
 * 
 * @param  M       Matrix
 * @return         C string of resolutions
 */
template <class T> inline static const char* 
ResToCString (const Matrix<T>& M) {
	return ResToString(M).c_str();
}

template<class T> inline std::ostream&
operator<< (std::ostream& os, const Matrix<T>& v) {
    size_t i,j;
    for (i = 0; i < v.Dim(0); ++i) {
        for (j = 0; j < v.Dim(1); ++j)
            TypeTraits<T>::print(os, v(i,j)) << " ";
        os << std::endl;
    }
    return os;
}


#endif
