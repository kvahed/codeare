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
//#include "Trigonometry.hpp"
//#include "Print.hpp"

/**
 * @brief    Absolute values 
 * 
 * @param  m Input
 * @return   Absolute values 
 */
template<class T> inline static Matrix<typename TypeTraits<T>::RT>
abs (const Matrix<T>& m) {
	Matrix<typename TypeTraits<T>::RT> res (size(m));
	for (size_t i = 0; i < numel(m); i++)
		res[i] = TypeTraits<T>::Abs(m[i]);
	return res;
}
 


/**
 * @brief    Arguments
 * 
 * @param  m Input
 * @return   Arguments
 */
template<class T> inline static Matrix<typename TypeTraits<T>::RT>
arg (const Matrix<T>& m) {
	Matrix<typename TypeTraits<T>::RT> res (size(m));
	for (size_t i = 0; i < numel(m); i++)
		res[i] = TypeTraits<T>::Arg(m[i]);
	return res;
}
 

/**
 * @brief    Real part
 * 
 * @param  m Input
 * @return   Real part
 */
template<class T> inline static Matrix<typename TypeTraits<T>::RT>
real (const Matrix<T>& m) {
	Matrix<typename TypeTraits<T>::RT> res (size(m));
	for (size_t i = 0; i < numel(m); i++)
		res[i] = TypeTraits<T>::Real(m[i]);
	return res;
}


/**
 * @brief    Imaginary part
 * 
 * @param  m Input
 * @return   Imaginary part
 */
template<class T> inline static Matrix<typename TypeTraits<T>::RT>
imag (const Matrix<T>& m) {
	Matrix<typename TypeTraits<T>::RT> res (size(m));
	for (size_t i = 0; i < numel(m); i++)
		res[i] = TypeTraits<T>::Imag(m[i]);
	return res;
}
 
template<class T> inline static Matrix<T>
conj (const MatrixType<T>& m) {
	Matrix<T> res (size(m));
	for (size_t i = 0; i < numel(m); i++)
		res[i] = TypeTraits<T>::Conj(m[i]);
	return res;
}

template<class T> inline static Matrix<T> conj (const Matrix<T>& m) {
	return (!m);
}

template <class T, class S> inline static Matrix<std::complex<T> > 
complex (const Matrix<T>& re, const Matrix<S>& im) {
    assert (numel(re) == numel(im));
    Matrix<std::complex<T> > ret (size(re));
    for (size_t i = 0; i < numel(re); ++i)
        ret[i] = std::complex<T>(re[i],im[i]);
    return ret;
}

template <class T, class S> inline static Matrix<std::complex<T> > 
	complex2 (const Matrix<T>& mag, const Matrix<S>& arg) {
    assert (numel(mag) == numel(arg));
    Matrix<std::complex<T> > ret (size(arg));
    for (size_t i = 0; i < numel(arg); ++i)
        ret[i] = std::polar(mag[i],arg[i]);
    return ret;
}
template <class T, class S> inline static Matrix<std::complex<T> > 
complex2 (const T mag, const Matrix<S>& arg) {
    Matrix<std::complex<T> > ret (size(arg));
    for (size_t i = 0; i < numel(arg); ++i)
        ret[i] = std::polar(mag[i],arg[i]);
    return ret;
}


#endif

