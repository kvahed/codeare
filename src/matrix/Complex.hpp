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

#ifndef __COMPLEX_HPP__
#define __COMPLEX_HPP__

#include "TypeTraits.hpp"

#include <stdlib.h>
#include <cstdlib>

template<class T>
static inline bool fp_type (const T t) {
	return (typeid(T) == typeid(float) || typeid(T) == typeid(double) ||
			typeid(T) ==  typeid(cxfl) || typeid(T) ==   typeid(cxdb) );
}

template<class T>
static inline bool is_complex (const T t) {
	return (typeid(T) == typeid(cxfl) || typeid(T) == typeid(cxdb) );
}

template<class T>
static inline bool is_singlep (const T t) {
	return (typeid(T) == typeid(cxfl) || typeid(T) == typeid(float) );
}

template<class T>
static inline bool is_doublep (const T t) {
	return (typeid(T) == typeid(cxdb) || typeid(T) == typeid(double) );
}

template<class T>
static inline bool i_type (const T t) {
	return (typeid(T) == typeid(short) || typeid(T) == typeid(long) ||
			typeid(T) ==  typeid(bool) || typeid(T) ==   typeid(int) );
}

template<class T>
static inline bool is_unsigned (const T t) {
	return (typeid(T) == typeid(unsigned int) || typeid(T) == typeid(unsigned short) ||
			typeid(T) == typeid(unsigned long) || typeid(T) == typeid(size_t));
}


template<class T>
struct CompTraits;

template<>
struct CompTraits<float> {

	typedef float type;

	inline static bool less_or_equal    (const type& a, const type& b) { return a <= b; }
	inline static bool less             (const type& a, const type& b) { return a <  b; }
	inline static bool greater_or_equal (const type& a, const type& b) { return a >= b; }
	inline static bool greater          (const type& a, const type& b) { return a >  b; }
	inline static bool logical_or       (const type& a, const type& b) { return a || b; }
	inline static bool logical_and      (const type& a, const type& b) { return a && b; }

};

template<>
struct CompTraits<double> {

	typedef double type;

	inline static bool less_or_equal    (const type& a, const type& b) { return a <= b; }
	inline static bool less             (const type& a, const type& b) { return a <  b; }
	inline static bool greater_or_equal (const type& a, const type& b) { return a >= b; }
	inline static bool greater          (const type& a, const type& b) { return a >  b; }
	inline static bool logical_or       (const type& a, const type& b) { return a || b; }
	inline static bool logical_and      (const type& a, const type& b) { return a && b; }

};

template<>
struct CompTraits<cxfl> {

	typedef cxfl type;

	inline static bool less_or_equal    (const type& a, const type& b) { return std::abs(a) <= std::abs(b); }
	inline static bool less             (const type& a, const type& b) { return std::abs(a) <  std::abs(b); }
	inline static bool greater_or_equal (const type& a, const type& b) { return std::abs(a) >= std::abs(b); }
	inline static bool greater          (const type& a, const type& b) { return std::abs(a) >  std::abs(b); }
	inline static bool logical_or       (const type& a, const type& b) { return std::abs(a) || std::abs(b); }
	inline static bool logical_and      (const type& a, const type& b) { return std::abs(a) && std::abs(b); }

};

template<>
struct CompTraits<cxdb> {

	typedef cxdb type;

	inline static bool less_or_equal    (const type& a, const type& b) { return std::abs(a) <= std::abs(b); }
	inline static bool less             (const type& a, const type& b) { return std::abs(a) <  std::abs(b); }
	inline static bool greater_or_equal (const type& a, const type& b) { return std::abs(a) >= std::abs(b); }
	inline static bool greater          (const type& a, const type& b) { return std::abs(a) >  std::abs(b); }
	inline static bool logical_or       (const type& a, const type& b) { return std::abs(a) || std::abs(b); }
	inline static bool logical_and      (const type& a, const type& b) { return std::abs(a) && std::abs(b); }

};

template<>
struct CompTraits<unsigned char> {

	typedef unsigned char type;

	inline static bool less_or_equal    (const type& a, const type& b) { return a <= b; }
	inline static bool less             (const type& a, const type& b) { return a <  b; }
	inline static bool greater_or_equal (const type& a, const type& b) { return a >= b; }
	inline static bool greater          (const type& a, const type& b) { return a >  b; }
	inline static bool logical_or       (const type& a, const type& b) { return a || b; }
	inline static bool logical_and      (const type& a, const type& b) { return a && b; }

};

template<>
struct CompTraits<size_t> {

	typedef size_t type;

	inline static bool less_or_equal    (const type& a, const type& b) { return a <= b; }
	inline static bool less             (const type& a, const type& b) { return a <  b; }
	inline static bool greater_or_equal (const type& a, const type& b) { return a >= b; }
	inline static bool greater          (const type& a, const type& b) { return a >  b; }
	inline static bool logical_or       (const type& a, const type& b) { return a || b; }
	inline static bool logical_and      (const type& a, const type& b) { return a && b; }

};


#endif
