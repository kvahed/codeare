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


#include <stdlib.h>
#include <complex>

typedef std::complex<float>  cxfl;
typedef std::complex<double> cxdb;

inline double cconj (double d) {return d;}
inline float  cconj (float  f) {return f;}
inline short  cconj (short  s) {return s;}
inline long   cconj (long   l) {return l;}
inline cxdb   cconj (cxdb  cd) {return std::conj(cd);}
inline cxfl   cconj (cxfl  cf) {return std::conj(cf);}

inline double creal (double d) {return d;}
inline float  creal (float  f) {return f;}
inline short  creal (short  s) {return s;}
inline long   creal (long   l) {return l;}
inline double creal (cxdb  cd) {return cd.real();}
inline float  creal (cxfl  cf) {return cf.real();}

inline double cimag (double d) {return 0.0;}
inline float  cimag (float  f) {return 0.0;}
inline short  cimag (short  s) {return 0;}
inline long   cimag (long   l) {return 0;}
inline double cimag (cxdb  cd) {return cd.imag();}
inline float  cimag (cxfl  cf) {return cf.imag();}

inline double cabs  (double d) {return fabs(d);}
inline float  cabs  (float  f) {return fabs(f);}
inline short  cabs  (short  s) {return fabs(s);}
inline long   cabs  (long   l) {return fabs(l);}
inline double cabs  (cxdb  cd) {return std::abs(cd);}
inline float  cabs  (cxfl  cf) {return std::abs(cf);}

inline double carg  (double d) {return 0.0;}
inline float  carg  (float  f) {return 0.0;}
inline short  carg  (short  s) {return 0;}
inline long   carg  (long   l) {return 0;}
inline double carg  (cxdb  cd) {return std::arg(cd);}
inline float  carg  (cxfl  cf) {return std::arg(cf);}

