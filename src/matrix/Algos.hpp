
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

#ifndef __ALGOS_HPP__
#define __ALGOS_HPP__

#include "IO.hpp"


/**
 * @brief    Number of non-zero elements
 *
 * @param M  Matrix
 * @return   Number of non-zero elements of matrix M
 */
template <class T> inline size_t 
nnz (const Matrix<T>& M) {
	
	size_t nz   = 0;
	T      zero = T(0);
	
	for (int i = 0; i < M.Size(); i++)
		if (M[i] != T(0))
			nz++;
	
	return nz;
	
}


/**
 * @brief    Is matrix 1D?
 *
 * @param M  Matrix
 * @return   1D?
 */
template <class T>  inline bool 
Is1D (const Matrix<T>& M) {
	
	return IsXD(M, 1);
	
}


/**
 * @brief    Is matrix 2D?
 *
 * @param M  Matrix
 * @return   2D?
 */
template <class T>  inline bool 
Is2D (const Matrix<T>& M) {
	
	return IsXD(M, 2);
	
}


/**
 * @brief    Is matrix 3D?
 *
 * @param M  Matrix
 * @return   3D?
 */
template <class T>  inline bool 
Is3D (const Matrix<T>& M) {
	
	return IsXD(M, 3);
	
}



/**
 * @brief    Is matrix 4D?
 *
 * @param M  Matrix
 * @return   4D?
 */
template <class T>  inline bool 
Is4D (const Matrix<T>& M) {
	
	return IsXD(M, 4);
	
}


/**
 * @brief     Is matrix X-dimensional?
 *
 * @param  M  Matrix
 * @param  d  Dimension
 * @return    X-dimensional?
 */
template <class T>  inline bool 
IsXD (const Matrix<T>& M, const size_t d) {
	
	size_t l = 0;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		if (M.Dim(i) > 1) l++;
	
	return (l == d);
	
}


/**
 * @brief       All elements zero?
 * 
 * @param  M    Matrix
 * @return      All elements zero?
 */
template <class T>  inline bool 
IsZero (const Matrix<T>& M) {
	
	for (size_t i = 0; i < M.Size(); i++)
		if (M[i] != T(0)) return false;
	
	return true;
	
}


/**
 * @brief       Empty matrix?
 * 
 * @param  M    Matrix
 * @return      Empty?
 */
template <class T>  inline bool
IsEmpty (const Matrix<T>& M) {
	
	return (numel(M) == 1);
	
}


/**
 * @brief       Sum of squares over a dimension
 * 
 * @param  M    Matrix
 * @param  d    Dimension
 * @return      Sum of squares
 */
template <class T> inline  Matrix<T> 
SOS (const Matrix<T>& M, const size_t d) {
	
	assert (M.Dim(d) > 1);
	
	Matrix<T> res = M ^ 2;
	return sum (res, d);
  
}


/**
 * @brief          Get rid of unused dimension 
 *
 * @param  M       Matrix
 * @return         Squeezed matrix
 */
template <class T> inline  Matrix<T>
squeeze (Matrix<T>& M) {
	
	size_t found = 0;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		if (M.Dim(i) > 1) {
			M.Res(found)   = M.Res(i);
			M.Dim(found++) = M.Dim(i);
		}
	
	for (size_t i = found; i < INVALID_DIM; i++) {
		M.Dim(i) = 1;
		M.Res(i) = 1.0;
	}

	return M;
	
}


/**
 * @brief     Sum along a dimension
 *
 * @param  M  Matrix
 * @param  d  Dimension
 * @return    Sum of M along dimension d
 */
template <class T> inline  Matrix<T>
sum (Matrix<T>& M, const size_t d) {
	
	Matrix<size_t> sz = size(M);
	size_t dim = M.Dim(d);
	Matrix<T> res;

	assert (d < INVALID_DIM);
	
	// No meaningful sum over particular dimension
	if (M.Dim(d) == 1) return res;
	
	// Empty? allocation 
	if (IsEmpty(M))    return res;
	
	// Inner size 
	size_t insize = 1;
	for (size_t i = 0; i < d; i++)
		insize *= M.Dim(i);
	
	// Outer size
	size_t outsize = 1;
	for (size_t i = d+1; i < INVALID_DIM; i++)
		outsize *= M.Dim(i);
	

	// Adjust size vector and allocate
	sz [d] = 1;

	res = Matrix<T>((size_t*)&sz[0]);

	// Sum
#pragma omp parallel default (shared) 
	{
		
		size_t tid      = omp_get_thread_num();
		
		for (size_t i = 0; i < outsize; i++) {
			
#pragma omp for
			
			for (size_t j = 0; j < insize; j++) {
				res[i*insize + j] = T(0);
				for (size_t k = 0; k < dim; k++)
					res[i*insize + j] += M[i*insize*dim + j + k*insize];
			}
			
		}
			
	}
	
	return res;
	
}



/**
 * @brief     Highest dimension unequal 1
 * 
 * @param  M  Matrix
 * @return    Highest non-one dimension
 */
template <class T> inline  size_t
ndims (const Matrix<T>& M) {
	
	size_t nd = 0;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		nd  = (M.Dim(i) > 1) ? i : nd;
	
	return nd;
	
}


/**
 * @brief           RAM occupied
 *
 * @param  M        Matrix
 * @return          Size in RAM in bytes.
 */
template <class T>  inline size_t
SizeInRAM          (const Matrix<T>& M) {
	
	return numel(M) * sizeof (T);
	
}


/**
 * @brief           Get the number of matrix cells, i.e. dim_0*...*dim_16.
 *
 * @param   M       Matrix
 * @return          Number of cells.
 */
template <class T>  size_t
numel               (const Matrix<T>& M) {
	
	size_t s = 1;
    
	for (size_t i = 0; i < INVALID_DIM; i++)
		s *= M.Dim(i);
    
	return s;
	
}


/**
 * @brief           Get vector of dimensions
 *
 * @param   M       Matrix
 * @return          Number of cells.
 */
template <class T>  Matrix<size_t>
size               (const Matrix<T>& M) {
	
	Matrix<size_t> res (1,INVALID_DIM);
	size_t ones = 0;
    
	for (size_t i = 0; i < INVALID_DIM; i++) {
		
		res[i] = M.Dim(i);
		ones   = (res[i] == 1) ? ones + 1 : ones = 0;
		
	}
	
	res.Resize (1,INVALID_DIM-ones);

	return res;
	
}


/**
 * @brief           Get size of a dimension
 *
 * @param   M       Matrix
 * @param   d       Dimension
 * @return          Number of cells.
 */
template <class T>  size_t
size               (const Matrix<T>& M, const size_t& d) {
	
	return M.Dim(d);
	
}



/**
 * @brief           Get length
 *
 * @param   M       Matrix
 * @return          Length
 */
template <class T>  size_t
length             (const Matrix<T>& M) {
	
	size_t l = 1;
	size_t cd;

	for (size_t i = 0; i < INVALID_DIM; i++) 
		l = (l > cd = M.Dim(i)) ? l : cd;

	return M.Dim(0);
	
}



/**
 * @brief           Get Width
 *
 * @param   M       Matrix
 * @return          Width
 */
template <class T>  size_t
width             (const Matrix<T>& M) {
	
	return M.Dim(1);
	
}



/**
 * @brief           Get height
 *
 * @param   M       Matrix
 * @return          Height
 */
template <class T>  size_t
height             (const Matrix<T>& M) {
	
	return M.Dim(0);
	
}



/**
 * @brief           Maximal element
 *
 * @param  M        Matrix
 * @return          Maximum
 */
template <class T> static inline T
max (const Matrix<T>& M) {

	T max = M[0];

	for (size_t i = 0; i < numel(M); i++)
		if (M[i] > max)
			max = M[i];

	return max;
	
}


/**
 * @brief           Minimal element
 *
 * @param  M        Matrix
 * @return          Minimum
 */
template <class T> static inline T
min (const Matrix<T>& M) {

	T max = M[0];

	for (size_t i = 0; i < numel(M); i++)
		if (M[i] < max)
			max = M[i];

	return max;
	
}


/**
 * @brief          FLip up down
 * 
 * @param          Matrix
 * @return         Flipped matrix
 */
/*
template <class T> inline static Matrix<T>
flipud (const Matrix<T>& M)  {

	size_t scol = size(M,0);
	size_t ncol = numel (M)/scol;

	Matrix<T> res (M.Dim());

	for (size_t i = 0; i < ncol; i++)
		reverse_copy (myints, myints+9, myvector.begin());

	return res;

}
*/
#endif 
