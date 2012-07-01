
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

#ifndef __ALGOS_HPP__
#define __ALGOS_HPP__

#include "IO.hpp"
#include "math.h"
#include "limits.h"

template <bool> struct static_assert;
template <> struct static_assert<true> { };

template<class T> inline bool 
is_nan (T const& x) {
    static_cast<void>(sizeof(static_assert<std::numeric_limits<T>::has_quiet_NaN>));
    return std::numeric_limits<T>::has_quiet_NaN and (x != x);
}

template <class T> inline bool 
is_inf (T const& x) {
    static_cast<void>(sizeof(static_assert<std::numeric_limits<T>::has_infinity>));
    return x == std::numeric_limits<T>::infinity() or x == -std::numeric_limits<T>::infinity();
}


/**
 * @brief    Number of non-zero elements
 *
 * Usage:
 * @code
 *   Matrix<cxfl> m = rand<double> (8,4,9,1,4);
 * 
 *   size_t nonz = nnz (M); // 1152
 * @endcode
 *
 * @param M  Matrix
 * @return   Number of non-zero elements of matrix M
 */
template <class T> inline static size_t 
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
template <class T>  inline static bool 
Is1D (const Matrix<T>& M) {
	
	return IsXD(M, 1);
	
}


/**
 * @brief    Is matrix 2D?
 *
 * @param M  Matrix
 * @return   2D?
 */
template <class T>  inline static bool 
Is2D (const Matrix<T>& M) {
	
	return IsXD(M, 2);
	
}


/**
 * @brief    Is 2D square matrix?
 *
 * @param M  Matrix
 * @return   2D?
 */
template <class T>  inline static bool 
IsSquare (const Matrix<T>& M) {
	
	return Is2D(M) && size(M,0) == size(M,1);
	
}


/**
 * @brief    Is matrix 3D?
 *
 * @param M  Matrix
 * @return   3D?
 */
template <class T>  inline static bool 
Is3D (const Matrix<T>& M) {
	
	return IsXD(M, 3);
	
}



/**
 * @brief    Is matrix 4D?
 *
 * @param M  Matrix
 * @return   4D?
 */
template <class T>  inline static bool 
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
template <class T>  inline static bool 
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
template <class T>  inline static bool 
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
template <class T> inline static bool
IsEmpty (const Matrix<T>& M) {
	
	return (numel(M) == 1);
	
}


/**
 * @brief       Which elements are NaN
 *
 * @param  M    Matrix
 * @return      Matrix of booleans true where NaN
 */
template <class T> inline static Matrix<bool> 
isnan (const Matrix<T>& M) {

    Matrix<bool> res (M.Dim());
	size_t i = numel(M);

	while (i--)
		res[i] = (std::isinf(creal(M[i]))||std::isinf(cimag(M[i])));
	
    return res;

}


/**
 * @brief       Which elements are Inf
 *
 * @param  M    Matrix
 * @return      Matrix of booleans true where inf
 */
template <class T> inline static Matrix<bool> 
isinf (const Matrix<T>& M) {

    Matrix<bool> res (M.Dim());
	size_t i = numel(M);

	while (i--)
		res.Dat().at(i) = (std::isinf(creal(M.At(i)))||std::isinf(cimag(M.At(i))));
	
    return res;

}


/**
 * @brief       Which elements are Inf
 *
 * @param  M    Matrix
 * @return      Matrix of booleans true where inf
 */
template <class T> inline static Matrix<bool> 
isfinite (const Matrix<T>& M) {

    Matrix<bool> res (M.Dim());
	size_t i = numel(M);

	while (i--)
		res.Dat()[i] = true;//is_nan(creal(M[i]));
			//()(!((bool)isnan(creal(M[i]))) && !isnan(cimag(M[i])) && !std::isinf(creal(M[i])) && !std::isinf(cimag(M[i])));
	
    return res;

}


/**
 * @brief       Make non finite elements (default T(0) else specify)
 *
 * @param  M    Matrix
 * @param  v    Optional value to which NaN and Inf elements are set. (default: T(0))
 * @return      Matrix stripped of NaN and Inf elements 
 */
template <class T> inline static Matrix<T> 
dofinite (const Matrix<T>& M, const T& v = 0) {

    Matrix<T> res (M.Dim());
	size_t i = numel(M);

	while (i--)
		res[i] = is_nan(creal(M[i])) ? v : M[i];
	//(std::isnan(ceal(M[i])) || std::isnan(cimag(M[i])) || std::isinf(creal(M[i])) || std::isinf(cimag(M[i]))) ? v : M[i];
	
    return res;

}


/**
 * @brief       Sum of squares over a dimension
 * 
 * Usage:
 * @code
 *   Matrix<cxfl> m   = rand<double> (8,7,6);
 *   m = sos (M,1); // dims (8,6);
 * @endcode
 *
 * @param  M    Matrix
 * @param  d    Dimension
 * @return      Sum of squares
 */
template <class T> inline static  Matrix<T> 
SOS (const Matrix<T>& M, const size_t d) {
	
	assert (M.Dim(d) > 1);
	
	Matrix<T> res = M ^ 2;
	return sum (res, d);
  
}


/**
 * @brief          Get rid of unused dimensions
 *
 * Usage:
 * @code
 *   Matrix<cxfl> m   = rand<double> (1,8,7,1,6);
 *   m = squeeze (m); // dims: (8,7,6); 
 * @endcode
 *
 * @param  M       Matrix
 * @return         Squeezed matrix
 */
template <class T> inline static Matrix<T>
squeeze (const Matrix<T>& M) {
	
	Matrix<T> res = M;
	size_t found = 0;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		if (size(res, i) > 1) {
			res.Res(found)   = M.Res(i);
			res.Dim(found++) = M.Dim(i);
		}
	
	for (size_t i = found; i < INVALID_DIM; i++) {
		res.Dim(i) = 1;
		res.Res(i) = 1.0;
	}
	
	return res;
	
}


/**
 * @brief     Sum along a dimension
 *
 * Usage:
 * @code
 *   Matrix<cxfl> m   = rand<double> (8,7,6);
 *   m = sum (m,0); // dims (7,6);
 * @endcode
 *
 * @param  M  Matrix
 * @param  d  Dimension
 * @return    Sum of M along dimension d
 */
template <class T> inline static Matrix<T>
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
 * Usage:
 * @code
 *   Matrix<cxfl> m   = rand<double> (8,7,6);
 *   size_t  nd       = ndims(m); // 2
 * @endcode
 *
 * @param  M  Matrix
 * @return    Highest non-one dimension
 */
template <class T> inline static size_t
ndims (const Matrix<T>& M) {
	
	size_t nd = 0;
	
	for (size_t i = 1; i < INVALID_DIM; i++)
		if (size(M,i) > 1)
			nd = i;
	
	return (nd + 1);
	
}





/**
 * @brief     Diagonal of biggest square matrix from top left
 * 
 * Usage:
 * @code
 *   Matrix<cxfl> m   = rand<double> (2,3);
 *   Matrix<cxfl> d   = diag (m);
 * @endcode
 *
 * @param  M  Matrix
 * @return    Highest non-one dimension
 */
template <class T> inline static Matrix<T>
diag (const Matrix<T>& M) {
	
	assert (Is2D(M));

	size_t sz = MIN(size(M,0),size(M,1));

	Matrix<T> res (sz,1);

	for (size_t i = 0; i < sz; i++)
		res(i) = M(i,i); 

	return res;
	
}


/**
 * @brief           RAM occupied
 *
 * @param  M        Matrix
 * @return          Size in RAM in bytes.
 */
template <class T> inline static size_t
SizeInRAM          (const Matrix<T>& M) {
	
	return numel(M) * sizeof (T);
	
}


/**
 * @brief           Get the number of matrix cells, i.e. dim_0*...*dim_16.
 *
 * @param   M       Matrix
 * @return          Number of cells.
 */
template <class T> inline static size_t
numel               (const Matrix<T>& M) {
	
	size_t s = 1;
    
	for (size_t i = 0; i < INVALID_DIM; i++)
		s *= M.Dim(i);
    
	return s;
	
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
template <class T>  inline static size_t
length             (const Matrix<T>& M) {
	
	size_t l = 1;

	for (size_t i = 0; i < INVALID_DIM; i++) 
		l = (l > size(M,i)) ? l : size(M,i);

	return l;
	
}



/**
 * @brief           Get Width
 *
 * @param   M       Matrix
 * @return          Width
 */
template <class T> inline static size_t
width             (const Matrix<T>& M) {
	
	return size(M,1);
	
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
template <class T> inline static T
max (const Matrix<T>& M) {

	T max = cabs(M[0]);

	for (size_t i = 0; i < numel(M); i++)
		if (cabs(M[i]) > max)
			max = M[i];

	return max;
	
}


/**
 * @brief           Minimal element
 *
 * @param  M        Matrix
 * @return          Minimum
 */
template <class T> inline static T
min (const Matrix<T>& M) {

	T min = cabs(M[0]);

	for (size_t i = 0; i < numel(M); i++)
		if (cabs(M[i]) < min)
			min = M[i];

	return min;
	
}



/**
 * @brief           Transpose
 *
 * @param  M        2D Matrix
 * @param  conj     Conjugate while transosing
 *
 * @return          Non conjugate transpose
 */
template <class T> inline static Matrix<T>
transpose (const Matrix<T>& M, const bool& conj = false) {

	assert (Is2D(M));
	size_t m, n, i, j;

	m = size(M,0);
	n = size(M,1);

	Matrix<T> res (n,m);

	for (j = 0; j < n; j++)
		for (i = 0; i < m; i++)
			res (j,i) = M(i,j);
	
	return (conj) ? conj(res) : res;

}


/**
 * @brief           Complex conjugate transpose
 *
 * @param  M        2D Matrix
 * @return          Complex conjugate transpose
 */
template <class T> inline static Matrix<T>
ctranspose (const Matrix<T>& M) {

	return transpose (M, true);

}


/**
 * @brief           MATLAB-like permute
 * 
 * @param   M       Input matrix 
 * @param   perm    New permuted dimensions
 * @return          Permuted matrix
 */
#include "Creators.hpp"
template <class T> inline static Matrix<T>
permute (const Matrix<T>& M, const Matrix<size_t>& perm) {
	
	// Check that perm only includes one number between 0 and INVALID_DIM once
	size_t ndnew = numel(perm), i, j;
	size_t ndold = ndims (M); 

	// Must have same number of dimensions
	assert (ndnew == ndold);

	// Every number between 0 and ndnew must appear exactly once
	vector<bool> occupied;
	occupied.resize(ndnew);
	for (i = 0; i < ndnew; i++) {
		assert (!occupied[perm[i]]);
		occupied [perm[i]] = true;
	}			

	// Old and new sizes
	Matrix<size_t> so = size (M);
	Matrix<size_t> sn = ones<size_t> (16,1);

	for (i = 0; i < ndnew; i++)
		sn[i] = so[perm[i]];
	
	// Allocate new matrix with permuted dimensions
	Matrix<T> res = zeros<T>(sn);

	// Relation of old to new indices
	size_t  d[16];
	size_t od[16];
	for (i = 0; i < ndnew; i++) od[i] = perm[i];
	for (     ; i <    16; i++)	od[i] =      i;
	
	// Copy data accordingly
	for (d[15] = 0; d[15] < size(res,15); d[15]++)
		for (d[14] = 0; d[14] < size(res,14); d[14]++)
			for (d[13] = 0; d[13] < size(res,13); d[13]++)
				for (d[12] = 0; d[12] < size(res,12); d[12]++)
					for (d[11] = 0; d[11] < size(res,11); d[11]++)
						for (d[10] = 0; d[10] < size(res,10); d[10]++)
							for (d[ 9] = 0; d[ 9] < size(res, 9); d[ 9]++)
								for (d[ 8] = 0; d[ 8] < size(res, 8); d[ 8]++)
									for (d[ 7] = 0; d[ 7] < size(res, 7); d[ 7]++)
										for (d[ 6] = 0; d[ 6] < size(res, 6); d[ 6]++)
											for (d[ 5] = 0; d[ 5] < size(res, 5); d[ 5]++)
												for (d[ 4] = 0; d[ 4] < size(res, 4); d[ 4]++)
													for (d[ 3] = 0; d[ 3] < size(res, 3); d[ 3]++)
														for (d[ 2] = 0; d[ 2] < size(res, 2); d[ 2]++)
															for (d[ 1] = 0; d[ 1] < size(res, 1); d[ 1]++)
																for (d[ 0] = 0; d[ 0] < size(res, 0); d[ 0]++) 
																	res (d[ 0],d[ 1],d[ 2],d[ 3],d[ 4],d[ 5],
																		 d[ 6],d[ 7],d[ 8],d[ 9],d[10],d[11],
																		 d[12],d[13],d[14],d[15]) =
																	  M (d[od[ 0]],d[od[ 1]],d[od[ 2]],d[od[ 3]], 
																		 d[od[ 4]],d[od[ 5]],d[od[ 6]],d[od[ 7]], 
																		 d[od[ 8]],d[od[ 9]],d[od[10]],d[od[11]],
																		 d[od[12]],d[od[13]],d[od[14]],d[od[15]]);

	return res;

																		
}


/**
 * @brief           Create new vector
 *                  and copy the data into the new vector. If the target
 *                  is bigger, the remaining space is set 0. If it is 
 *                  smaller data is truncted.
 * 
 * @param   M       The matrix to resize
 * @param   sz      New length
 * @return          Resized vector
 */
template <class T> inline static Matrix<T>
resize (const Matrix<T>& M, const size_t& sz) {

	Matrix<T> res = zeros<T> (sz,1);
	size_t copysz = MIN(numel(M), sz) * sizeof (T);

	memcpy (&res[0], &(M[0]), copysz);

	return res;
	
}



/**
 * @brief           Create new matrix with the new dimensions 
 *                  and copy the data into the new matrix. If the target
 *                  is bigger, the remaining space is set 0. If it is 
 *                  smaller data is truncted.
 * 
 * @param   M       The matrix to resize
 * @param   sz      New dimension vector
 * @return          Resized copy
 */
template <class T> inline static Matrix<T>
resize (const Matrix<T>& M, Matrix<size_t> sz) {

	Matrix<T> res = zeros<T> (sz);
	size_t copysz = MIN(numel(M), numel(res)) * sizeof (T);

	typename vector<T>::iterator start = M.Dat().begin();
	typename vector<T>::iterator end   = M.Dat().begin() + copysz;
	
	res.Dat().assign (start, end);

	return res;
	
}


/**
 * @brief           Create new vector
 *                  and copy the data into the new vector. If the target
 *                  is bigger, the remaining space is set 0. If it is 
 *                  smaller data is truncted.
 * 
 * @param   M       The matrix to resize
 * @param   sc      New height
 * @param   sl      New width
 * @return          Resized vector
 */
template <class T> inline static Matrix<T>
resize (const Matrix<T>& M, const size_t& sc, const size_t& sl) {

	Matrix<size_t> sz (2,1);
	sz[0] = sc; sz[1] = sl;

	return resize (M, sz);
	
}


/**
 * @brief           Create new vector
 *                  and copy the data into the new vector. If the target
 *                  is bigger, the remaining space is set 0. If it is 
 *                  smaller data is truncted.
 * 
 * @param   M       The matrix to resize
 * @param   sc      New height
 * @param   sl      New width
 * @param   sc      New # of slices
 * @return          Resized vector
 */
template <class T> inline static Matrix<T>
resize (const Matrix<T>& M, const size_t& sc, const size_t& sl, const size_t& ss) {

	Matrix<size_t> sz (3,1);
	sz[0] = sc; sz[1] = sl; sz[2] = sc;

	return resize (M, sz);
	
}


/**
 * @brief           Get vector of dimensions
 *
 * @param   M       Matrix
 * @return          Dimension vector.
 */
template <class T>  inline static Matrix<size_t>
size               (const Matrix<T>& M) {
	
	Matrix<size_t> res (1,INVALID_DIM);
	size_t ones = 0;
    
	for (size_t i = 0; i < INVALID_DIM; i++) {
		
		res[i] = size(M, i);

		if (res[i] == 1)
			ones++;
		else
			ones = 0;
		
	}
	
	return resize (res, 1, INVALID_DIM-ones);
	
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
