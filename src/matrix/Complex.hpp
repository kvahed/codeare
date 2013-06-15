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

#include <stdlib.h>
#include <complex>
#include <typeinfo>

typedef std::complex<float>  cxfl;
typedef std::complex<double> cxdb;


static const std::type_info& float_type (typeid(float));
static const std::type_info& double_type (typeid(double));
static const std::type_info& cxfl_type (typeid(cxfl));
static const std::type_info& cxdb_type (typeid(cxdb));

static const std::type_info& short_type (typeid(short));
static const std::type_info& long_type (typeid(long));
static const std::type_info& bool_type (typeid(bool));
static const std::type_info& int_type (typeid(int));
static const std::type_info& size_t_type (typeid(size_t));

template<class T>
static inline bool fp_type (const T t) {
	return (typeid(T) == float_type || typeid(T) == double_type ||
			typeid(T) ==  cxfl_type || typeid(T) ==   cxdb_type );
}

template<class T>
static inline bool is_complex (const T t) {
	return (typeid(T) == cxfl_type || typeid(T) == cxdb_type );
}

template<class T>
static inline bool is_singlep (const T t) {
	return (typeid(T) == cxfl_type || typeid(T) == float_type );
}

template<class T>
static inline bool is_doublep (const T t) {
	return (typeid(T) == cxdb_type || typeid(T) == double_type );
}

template<class T>
static inline bool i_type (const T t) {
return (typeid(T) == short_type || typeid(T) == long_type ||
		typeid(T) ==  bool_type || typeid(T) ==   int_type );
}

inline double cconj (double d) {return d;}
inline float  cconj (float  f) {return f;}
inline short  cconj (short  s) {return s;}
inline long   cconj (long   l) {return l;}
inline long   cconj (size_t s) {return s;}
inline cxdb   cconj (cxdb  cd) {return std::conj(cd);}
inline cxfl   cconj (cxfl  cf) {return std::conj(cf);}

inline double creal (double d) {return d;}
inline float  creal (float  f) {return f;}
inline short  creal (short  s) {return s;}
inline long   creal (long   l) {return l;}
inline long   creal (size_t s) {return s;}
inline double creal (cxdb  cd) {return cd.real();}
inline float  creal (cxfl  cf) {return cf.real();}

inline double cimag (double d) {return 0.0;}
inline float  cimag (float  f) {return 0.0;}
inline short  cimag (short  s) {return 0;}
inline long   cimag (long   l) {return 0;}
inline long   cimag (size_t s) {return 0;}
inline double cimag (cxdb  cd) {return cd.imag();}
inline float  cimag (cxfl  cf) {return cf.imag();}

inline double cabs  (double d) {return fabs(d);}
inline float  cabs  (float  f) {return fabs(f);}
inline short  cabs  (short  s) {return fabs(s);}
inline long   cabs  (long   l) {return fabs(l);}
inline long   cabs  (int    i) {return fabs(i);}
inline long   cabs  (size_t s) {return s;}
inline double cabs  (cxdb  cd) {return std::abs(cd);}
inline float  cabs  (cxfl  cf) {return std::abs(cf);}

inline double carg  (double d) {return 0.0;}
inline float  carg  (float  f) {return 0.0;}
inline short  carg  (short  s) {return 0;}
inline long   carg  (long   l) {return 0;}
inline long   carg  (size_t s) {return 0;}
inline double carg  (cxdb  cd) {return std::arg(cd);}
inline float  carg  (cxfl  cf) {return std::arg(cf);}

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

	inline static bool less_or_equal    (const type& a, const type& b) { return abs(a) <= abs(b); }
	inline static bool less             (const type& a, const type& b) { return abs(a) <  abs(b); }
	inline static bool greater_or_equal (const type& a, const type& b) { return abs(a) >= abs(b); }
	inline static bool greater          (const type& a, const type& b) { return abs(a) >  abs(b); }
	inline static bool logical_or       (const type& a, const type& b) { return abs(a) || abs(b); }
	inline static bool logical_and      (const type& a, const type& b) { return abs(a) && abs(b); }

};

template<>
struct CompTraits<cxdb> {

	typedef cxdb type;

	inline static bool less_or_equal    (const type& a, const type& b) { return abs(a) <= abs(b); }
	inline static bool less             (const type& a, const type& b) { return abs(a) <  abs(b); }
	inline static bool greater_or_equal (const type& a, const type& b) { return abs(a) >= abs(b); }
	inline static bool greater          (const type& a, const type& b) { return abs(a) >  abs(b); }
	inline static bool logical_or       (const type& a, const type& b) { return abs(a) || abs(b); }
	inline static bool logical_and      (const type& a, const type& b) { return abs(a) && abs(b); }

};

template<>
struct CompTraits<unsigned short> {

	typedef bool type;

	inline static bool less_or_equal    (const type& a, const type& b) { return a <= b; }
	inline static bool less             (const type& a, const type& b) { return a <  b; }
	inline static bool greater_or_equal (const type& a, const type& b) { return a >= b; }
	inline static bool greater          (const type& a, const type& b) { return a >  b; }
	inline static bool logical_or       (const type& a, const type& b) { return a || b; }
	inline static bool logical_and      (const type& a, const type& b) { return a && b; }

};



#endif
