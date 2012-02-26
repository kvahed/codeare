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

#ifndef __CX_HPP__
#define __CX_HPP__

#include "Matrix.hpp"


/**
 * @brief    Absolute values 
 * 
 * @param  m Input
 * @return   Absolute values 
 */
Matrix<double> 
abs (const Matrix<cxdb>& m) {

	Matrix<double> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
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
Matrix<float> 
abs (const Matrix<cxfl>& m) {

	Matrix<float> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
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
Matrix<float> 
abs (const Matrix<float>& m) {
	
	Matrix<float> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
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
Matrix<double> 
abs (const Matrix<double>& m) {
	
	Matrix<float> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
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
Matrix<double> 
arg (const Matrix<cxdb>& m) {

	Matrix<double> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
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
Matrix<float> 
arg (const Matrix<cxfl>& m) {

	Matrix<float> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
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
Matrix<float> 
arg (const Matrix<float>& m) {
	
	Matrix<float> res(m);
	return res;
	
}


/**
 * @brief    Arguments
 * 
 * @param  m Input
 * @return   Arguments
 */
Matrix<double> 
arg (const Matrix<double>& m) {
	
	Matrix<double> res(m);
	return res;
	
}


/**
 * @brief    Real part
 * 
 * @param  m Input
 * @return   Real part
 */
Matrix<double> 
real (const Matrix<cxdb>& m) {

	Matrix<double> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
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
Matrix<float> 
real (const Matrix<cxfl>& m) {

	Matrix<float> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
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
Matrix<float> 
real (const Matrix<float>& m) {
	
	Matrix<float> res(m);
	return res;
	
}


/**
 * @brief    Dummy real
 * 
 * @param  m Input
 * @return   Real part
 */
Matrix<double> 
real (const Matrix<double>& m) {
	
	Matrix<double> res(m);
	return res;
	
}



/**
 * @brief    Imaginary part
 * 
 * @param  m Input
 * @return   Imaginary part
 */
Matrix<double> 
imag (const Matrix<cxdb>& m) {

	Matrix<double> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
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
Matrix<float> 
imag (const Matrix<cxfl>& m) {

	Matrix<float> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
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
Matrix<float> 
imag (const Matrix<float>& m) {
	
	Matrix<float> res(m);
	return res;
	
}



/**
 * @brief    Dummy for double
 * 
 * @param  m Input
 * @return   Imaginary part
 */
Matrix<double> 
imag (const Matrix<double>& m) {
	
	Matrix<double> res(m);
	return res;
	
}


/**
 * @brief           Complex conjugate (no transposition)
 *
 * @param   m       Input 
 * @return          Complex conjugate 
 */
template<class T> Matrix<T>
conj (const Matrix<T>& m) {
	
	Matrix<T> res = m;
	
		if (typeid (T) == typeid (cxfl) || typeid (T) == typeid (cxdb)) {
			
#pragma omp parallel default (shared) 
			{
				
#pragma omp for
				
				for (size_t i = 0; i < m.Size(); i++)
					res[i] = cconj(m[i]);
				
			}		
			
		}
		
		return res;
		
}

#endif

