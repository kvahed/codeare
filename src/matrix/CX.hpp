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

#ifndef __CX_HPP__
#define __CX_HPP__

#include "Matrix.hpp"
#include "Complex.hpp"
#include "Creators.hpp"
#include "Algos.hpp"
#include "Trigonometry.hpp"

/**
 * @brief    Absolute values 
 * 
 * @param  m Input
 * @return   Absolute values 
 */
inline static Matrix<double> 
abs (const Matrix<cxdb>& m) {

	Matrix<double> res (size(m));

#pragma omp parallel default (shared) 
	{
#pragma omp for
		
		for (size_t i = 0; i < numel(m); i++)
			res[i] = cabs(m[i]);
		
	}		
	
	return res;
	
}
 

/**
 * @brief    Absolute values 
 * 
 * @param  m Input
 * @return   Absolute values 
 */
inline static Matrix<float> 
abs (const Matrix<cxfl>& m) {

	Matrix<float> res = zeros<float> (size(m));
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < numel(m); i++)
			res[i] = cabs(m[i]);
		
	}		
	
	return res;
	
}


/**
 * @brief    Absolute values 
 * 
 * @param  m Input
 * @return   Absolute values 
 */
inline static Matrix<float> 
abs (const Matrix<float>& m) {
	
	Matrix<float> res = zeros<float> (size(m));
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < numel(m); i++)
			res[i] = cabs(m[i]);
		
	}		
	
	return res;
	
}


/**
 * @brief    Absolute values 
 * 
 * @param  m Input
 * @return   Absolute values 
 */
inline static Matrix<double> 
abs (const Matrix<double>& m) {
	
	Matrix<double> res = zeros<double> (size(m));
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < numel(m); i++)
			res[i] = cabs(m[i]);
		
	}		
	
	return res;
	
}


/**
 * @brief    Arguments
 * 
 * @param  m Input
 * @return   Arguments
 */
inline static Matrix<double> 
arg (const Matrix<cxdb>& m) {

	Matrix<double> res = zeros<double> (size(m));
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < numel(m); i++)
			res[i] = carg(m[i]);
		
	}		
	
	return res;
	
}
 

/**
 * @brief    Arguments
 * 
 * @param  m Input
 * @return   Arguments
 */
inline static Matrix<float> 
arg (const Matrix<cxfl>& m) {

	Matrix<float> res = zeros<float> (size(m));
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < numel(m); i++)
			res[i] = carg(m[i]);
		
	}		
	
	return res;
	
}


/**
 * @brief    Arguments
 * 
 * @param  m Input
 * @return   Arguments
 */
inline static Matrix<float> 
arg (const Matrix<float>& m) {
	
	return 	zeros<float> (size(m));
	
}


/**
 * @brief    Arguments
 * 
 * @param  m Input
 * @return   Arguments
 */
inline static Matrix<double> 
arg (const Matrix<double>& m) {
	
	return zeros<double> (size(m));
	
}


/**
 * @brief    Real part
 * 
 * @param  m Input
 * @return   Real part
 */
inline static Matrix<double> 
real (const Matrix<cxdb>& m) {

	Matrix<double> res = zeros<double> (size(m));
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < numel(m); i++)
			res[i] = creal(m[i]);
		
	}		
	
	return res;
	
}
 

/**
 * @brief    Real part
 * 
 * @param  m Input
 * @return   Real part
 */
inline static Matrix<float> 
real (const Matrix<cxfl>& m) {

	Matrix<float> res = zeros<float> (size(m));
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < numel(m); i++)
			res[i] = creal(m[i]);
		
	}		
	
	return res;
	
}


/**
 * @brief    Dummy Real 
 * 
 * @param  m Input
 * @return   Real part
 */
inline static Matrix<float> 
real (const Matrix<float>& m) {
	
	return Matrix<float>(m);
	
}


/**
 * @brief    Dummy real
 * 
 * @param  m Input
 * @return   Real part
 */
inline static Matrix<double> 
real (const Matrix<double>& m) {
	
	return Matrix<double>(m);
	
}



/**
 * @brief    Imaginary part
 * 
 * @param  m Input
 * @return   Imaginary part
 */
inline static Matrix<double> 
imag (const Matrix<cxdb>& m) {

	Matrix<double> res = zeros<double> (size(m));
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < numel(m); i++)
			res[i] = cimag(m[i]);
		
	}		
	
	return res;
	
}
 

/**
 * @brief    Imaginary part
 * 
 * @param  m Input
 * @return   Imaginary part
 */
inline static Matrix<float> 
imag (const Matrix<cxfl>& m) {

	Matrix<float> res = zeros<float> (size(m));
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < numel(m); i++)
			res[i] = cimag(m[i]);
		
	}		
	
	return res;
	
}


/**
 * @brief    Dummy for float
 * 
 * @param  m Input
 * @return   Imaginary part
 */
inline static Matrix<float> 
imag (const Matrix<float>& m) {
	
	return 	zeros<double>(size(m));
	
}



/**
 * @brief    Dummy for double
 * 
 * @param  m Input
 * @return   Imaginary part
 */
inline static Matrix<double> 
imag (const Matrix<double>& m) {
	
	return zeros<double>(size(m));
	
}


/**
 * @brief           Complex conjugate (no transposition)
 *
 * @param   m       Input 
 * @return          Complex conjugate 
 */
template<class T> inline static Matrix<T>
conj (const Matrix<T>& m) {
	
	Matrix<T> res = m;
	
	if (typeid (T) == typeid (cxfl) || typeid (T) == typeid (cxdb)) {

#pragma omp parallel default (shared)
		{
			
#pragma omp for
			
			for (size_t i = 0; i < numel(m); i++)
				res[i] = cconj(res[i]);

		}
		
	}

	return res;

}


using namespace codeare::matrix::arithmetic;

template <class T, class S> inline static Matrix<std::complex<T> > 
complex (const Matrix<T>& re, const Matrix<S>& im) {
    assert (numel(re) == numel(im));
    Matrix<std::complex<T> > ret (vsize(re));
#pragma omp for
    for (size_t i = 0; i < numel(re); ++i)
        ret[i] = std::complex<T>(re[i],im[i]);
    return ret;
}

template <class T, class S> inline static Matrix<std::complex<T> > 
complex2 (const Matrix<T>& mag, const Matrix<S>& arg) {
    assert (numel(mag) == numel(arg));
    Matrix<std::complex<T> > ret (vsize(arg));
#pragma omp for
    for (size_t i = 0; i < numel(arg); ++i)
        ret[i] = std::complex<T>(mag[i]*cos(arg[i]),mag[i]*sin(arg[i]));
    return ret;
}
template <class T, class S> inline static Matrix<std::complex<T> > 
complex2 (const T mag, const Matrix<S>& arg) {
    Matrix<std::complex<T> > ret (vsize(arg));
#pragma omp for
    for (size_t i = 0; i < numel(arg); ++i)
        ret[i] = std::complex<T>(mag*cos(arg[i]),mag*sin(arg[i]));
    return ret;
}


#endif

