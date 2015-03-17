
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

#include "Matrix.hpp"
#include "math.h"
#include <limits>
#include <vector>
#include <algorithm>    // std::reverse


#if !defined(_MSC_VER) || _MSC_VER>1200
#include <boost/math/special_functions/fpclassify.hpp>

template<class T> inline static bool is_nan (T const& x) {
    return boost::math::isnan (x);
}
template <class T> inline static bool is_inf (T const& x) {
    return boost::math::isinf (x);
}
#endif



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
template <class T> inline static  size_t
nnz (const Matrix<T>& M) {
	
	size_t nz   = 0;
	
	for (size_t i = 0; i < M.Size(); ++i)
		if (M[i] != T(0))
			++nz;
	
	return nz;
	
}


/**
 * @brief     Is matrix X-dimensional?
 *
 * @param  M  Matrix
 * @param  d  Dimension
 * @return    X-dimensional?
 */
template <class T>  inline static  bool
isxd (const Matrix<T>& M, size_t d) {

	size_t l = 0;

	for (size_t i = 0; i < M.NDim(); ++i)
		if (M.Dim(i) > 1) ++l;

	return (l == d);

}

/**
 * @brief    Is matrix 1D?
 *
 * @param M  Matrix
 * @return   1D?
 */
template <class T>  inline static  bool
isvec (const Matrix<T>& M) {
	
	return isxd(M, 1);
	
}


/**
 * @brief    Is matrix 2D?
 *
 * @param M  Matrix
 * @return   2D?
 */
template <class T>  inline static  bool
is2d (const Matrix<T>& M) {
	
	return isxd(M, 2);
	
}


/**
 * @brief    Is 2D square matrix?
 *
 * @param M  Matrix
 * @return   2D?
 */
template <class T>  inline static  bool
issquare (const Matrix<T>& M) {
	
	return isxd(M, 2) && (size(M,0) == size(M,1));
	
}


/**
 * @brief    Is matrix 3D?
 *
 * @param M  Matrix
 * @return   3D?
 */
template <class T>  inline static  bool
is3d (const Matrix<T>& M) {
	
	return isxd(M, 3);
	
}



/**
 * @brief    Is matrix 4D?
 *
 * @param M  Matrix
 * @return   4D?
 */
template <class T>  inline static  bool
is4d (const Matrix<T>& M) {
	
	return isxd(M, 4);
	
}


/**
 * @brief       All elements zero?
 * 
 * @param  M    Matrix
 * @return      All elements zero?
 */
template <class T>  inline static  bool
iszero (const Matrix<T>& M) {
	
	for (size_t i = 0; i < M.Size(); ++i)
		if (M[i] != T(0)) return false;
	
	return true;
	
}


/**
 * @brief       Empty matrix?
 * 
 * @param  M    Matrix
 * @return      Empty?
 */
template <class T> inline static  bool
isempty (const Matrix<T>& M) {
	
	return (numel(M) == 1);
	
}


/**
 * @brief       Which elements are NaN
 *
 * @param  M    Matrix
 * @return      Matrix of booleans true where NaN
 */
template <class T> inline static  Matrix<cbool>
isinf (const Matrix<T>& M) {

    Matrix<cbool> res (M.Dim());
    for (size_t i = 0; i < res.Size(); ++i)
		res[i] = (is_inf(TypeTraits<T>::Real(M[i]))||is_inf(TypeTraits<T>::Imag(M[i])));
    return res;

}


/**
 * @brief       Which elements are Inf
 *
 * @param  M    Matrix
 * @return      Matrix of booleans true where inf
 */
template <class T> inline static  Matrix<cbool>
isnan (const Matrix<T>& M) {

    Matrix<cbool> res (M.Dim());
    for (size_t i = 0; i < res.Size(); ++i)
		res.Container()[i] = (is_nan(TypeTraits<T>::Imag(M[i]))||is_nan(TypeTraits<T>::Imag(M[i])));
    return res;

}


/**
 * @brief       Which elements are Inf
 *
 * @param  M    Matrix
 * @return      Matrix of booleans true where inf
 */
template <class T> inline static  Matrix<cbool>
isfinite (const Matrix<T>& M) {

    Matrix<cbool> res (M.Dim());
	size_t i = numel(M);

	
    return res;

}


/**
 * @brief       Make non finite elements (default T(0) else specify)
 *
 * @param  M    Matrix
 * @param  v    Optional value to which NaN and Inf elements are set. (default: T(0))
 * @return      Matrix stripped of NaN and Inf elements 
 */
#if !defined(_MSC_VER) || _MSC_VER>1200
template <class T> inline static  Matrix<T>
dofinite (const Matrix<T>& M, const T& v = 0) {

    Matrix<T> res (M.Dim());
	size_t i = numel(M);

	while (i--)
		res[i] = is_nan(TypeTraits<T>::Real(M[i])) ? v : M[i];
	
    return res;

}
#endif


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
template <class T> inline static  size_t
ndims (const Matrix<T>& M) {
	
	size_t nd = 0;
	
	for (size_t i = 1; i < M.NDim(); ++i)
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
template <class T> inline static  Matrix<T>
diag (const Matrix<T>& M) {
	
	assert (is2d(M));

	size_t sz = MIN(size(M,0),size(M,1));

	Matrix<T> res (sz,1);

	for (size_t i = 0; i < sz; ++i)
		res(i) = M(i,i); 

	return res;
	
}


/**
 * @brief           RAM occupied
 *
 * @param  M        Matrix
 * @return          Size in RAM in bytes.
 */
template <class T> inline static  size_t
SizeInRAM          (const Matrix<T>& M) {
	
	return numel(M) * sizeof (T);
	
}


/**
 * @brief           Get the number of matrix cells, i.e. dim_0*...*dim_16.
 *
 * @param   M       Matrix
 * @return          Number of cells.
 */
template <class T, paradigm P> inline static  size_t
numel               (const Matrix<T,P>& M) {
	return M.Size();
}


/**
 * @brief          All elements non-zero?
 *
 * @param   M      Matrix in question
 * @return         True if all elements non-zero
 */
template <class T> inline static bool
all (const Matrix<T>& M) {
	return (nnz(M) == numel(M));
}


/**
 * @brief           Get size of a dimension
 *
 * @param   M       Matrix
 * @param   d       Dimension
 * @return          Number of cells.
 */
template <class T>  size_t
size               (const Matrix<T>& M, size_t d) {
	return M.Dim(d);
}
template <class T>  size_t
size               (const Matrix<T,MPI>& M, size_t d) {
	return M.Dim(d);
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
template <class T> inline static  Matrix<T>
resize (const Matrix<T>& M, size_t sc, size_t sl) {

	Vector<size_t> sz (2);
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
 * @param   ss      New # of slices
 *
 * @return          Resized vector
 */
template <class T> inline static  Matrix<T>
resize (const Matrix<T>& M, size_t sc, size_t sl, size_t ss) {

	Matrix<size_t> sz (3,1);
	sz[0] = sc; sz[1] = sl; sz[2] = ss;

	return resize (M, sz);
	
}


/**
 * @brief           Get vector of dimensions
 *
 * @param   M       Matrix
 * @return          Dimension vector.
 */
template <class T,paradigm P>  inline static  Vector<size_t>
size               (const Matrix<T,P>& M) {
	
	Vector<size_t> ret (M.NDim());
	for (size_t i = 0; i < ret.size(); ++i)
		ret[i] = size(M,i);
	return ret;

}


/**
 * @brief           Get resolution of a dimension
 *
 * @param   M       Matrix
 * @param   d       Dimension
 * @return          Resolution
 */
template <class T>  size_t
resol               (const Matrix<T>& M, size_t d) {
	
	return M.Res(d);
	
}


/**
 * @brief           Get length
 *
 * @param   M       Matrix
 * @return          Length
 */
template <class T>  inline static  size_t
length             (const Matrix<T>& M) {
	
	size_t l = 1;

	for (size_t i = 0; i < M.NDim(); ++i)
		l = (l > size(M,i)) ? l : size(M,i);

	return l;
	
}



/**
 * @brief           Get Width
 *
 * @param   M       Matrix
 * @return          Width
 */
template <class T> inline static  size_t
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
 * @brief           Round down
 *
 * @param  M        Matrix
 * @return          Rounded down matrix
 */
template<class T> inline static Matrix<T>
floor (const Matrix<T>& M) {
	Matrix<T> res = M;
	for (size_t i = 0; i < numel(M); ++i)
		res[i] = floor ((float)res[i]);
	return res;
}


/**
 * @brief           Round up
 *
 * @param  M        Matrix
 * @return          Rounded up matrix
 */
template<class T> inline static Matrix<T>
ceil (const Matrix<T>& M) {
	Matrix<T> res = M;
	for (size_t i = 0; i < numel(M); ++i)
		res[i] = ceil (res[i]);
	return res;
}


/**
 * @brief           MATLAB-like round
 *
 * @param  M        Matrix
 * @return          Rounded matrix
 */
template<class T> inline static Matrix<T>
round (const Matrix<T>& M) {
	Matrix<T> res = M;
	for (size_t i = 0; i < numel(M); ++i)
		res[i] = ROUND (res[i]);
	return res;
}


/**
 * @brief           Maximal element
 *
 * @param  M        Matrix
 * @return          Maximum
 */
template <class T> inline static  T
m_max (const Matrix<T>& M) {

	T mx = std::abs(M[0]);

	for (size_t i = 0; i < numel(M); ++i)
		if (std::abs(M[i]) > mx)
			mx = M[i];

	return mx;
	
}

/**
 * @brief           Maximal element
 *
 * @param  M        Matrix
 * @return          Maximum
 */
template <class T> inline static  T
max (const Matrix<T>& M) {

	T mx = (T)INT_MIN;

	for (size_t i = 0; i < numel(M); ++i)
		if (M[i] > mx)
			mx = M[i];

	return mx;

}

/**
 * @brief           Maximal element
 *
 * @param  M        Matrix
 * @return          Maximum
 */
template <class T> inline static  T
min (const Matrix<T>& M) {

	T mn = (T)INT_MAX;

	for (size_t i = 0; i < numel(M); ++i)
		if (M[i] < mn)
			mn= M[i];

	return mn;

}

/**
 * @brief           Minimal element
 *
 * @param  M        Matrix
 * @return          Minimum
 */
template <class T> inline static  T
m_min (const Matrix<T>& M) {

	T mn = abs(M[0]);

	for (size_t i = 0; i < numel(M); ++i)
		if (abs(M[i]) < mn)
			mn = M[i];

	return mn;
	
}


/**
 * @brief           Transpose
 *
 * @param  M        2D Matrix
 * @param  c        Conjugate while transposing
 *
 * @return          Non conjugate transpose
 */
template <class T> inline static  Matrix<T>
transpose (const Matrix<T>& M, bool c = false) {

	assert (is2d(M));
	size_t m = size(M,0), n = size(M,1), i, j;
	Matrix<T> res (M);

	for (j = 0; j < n; ++j)
		for (i = 0; i < j; ++i)
			swapd(res(j,i),res(i,j));
	
	return c ? conj(res) : res;

}


/**
 * @brief           Complex conjugate transpose
 *
 * @param  M        2D Matrix
 * @return          Complex conjugate transpose
 */
template <class T> inline static  Matrix<T>
ctranspose (const Matrix<T>& M) {
	return transpose (M, true);
}


#include "Creators.hpp"

/*
 * @brief           Create new vector
 *                  and copy the data into the new vector. If the target
 *                  is bigger, the remaining space is set 0. If it is 
 *                  smaller data is truncted.
 * 
 * @param   M       The matrix to resize
 * @param   sz      New length
 * @return          Resized vector
 */
template <class T> inline static  Matrix<T>
resize (const Matrix<T>& M, size_t sz) {

	Matrix<T> res (sz,1);
	size_t copysz = MIN(numel(M), sz);

    typename Vector<T>::      iterator rb = res.Begin ();
    typename Vector<T>::const_iterator mb =   M.Begin ();
    
    std::copy (mb, mb+copysz, rb);

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
resize (const Matrix<T>& M, const Vector<size_t>& sz) {

	Matrix<T> res (sz);
	size_t copysz  = MIN(numel(M), numel(res));

    typename Vector<T>::      iterator rb = res.Begin ();
    typename Vector<T>::const_iterator mb =   M.Begin ();

    std::copy (mb, mb+copysz, rb);

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
sum (const Matrix<T>& M, size_t d) {
	
	Vector<size_t> sz = size(M);
	size_t        dim = sz[d];
	Matrix<T>     res;

	assert (d < M.NDim());
	
	// No meaningful sum over particular dimension
	if (dim == 1) 
		return res;
	
	// Empty? allocation 
	if (isempty(M))
		return res;
	
	// Inner size 
	size_t insize = 1;
	for (size_t i = 0; i < d; ++i)
		insize *= sz[i];
	
	// Outer size
	size_t outsize = 1;
	for (size_t i = d+1; i < MIN(M.NDim(),sz.size()); ++i)
		outsize *= sz[i];
	
	// Adjust size vector and allocate
	sz [d] = 1;
	res = Matrix<T>(sz);

	// Sum
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (int i = 0; i < outsize; ++i) {
			
			for (size_t j = 0; j < insize; ++j) {
				res[i*insize + j] = T(0);
				for (size_t k = 0; k < dim; ++k)
					res[i*insize + j] += M[i*insize*dim + j + k*insize];
			}
			
		}
		
	}

	return res;
	
}


/**
 * @brief     Sum of all elements
 *
 * Usage:
 * @code
 *   Matrix<cxfl> m   = rand<double> (8,7,6);
 *   m = sum (m);
 * @endcode
 *
 * @param  M  Matrix
 * @return    Sum of M along dimension d
 */
template <class T> inline static T
sum (const Matrix<T>& M) {
	return std::accumulate(M.Begin(),M.End(),(T)0);
}


/**
 * @brief     Product along a dimension
 *
 * Usage:
 * @code
 *   Matrix<cxfl> m   = rand<double> (8,7,6);
 *   m = prod (m,0); // dims (7,6);
 * @endcode
 *
 * @param  M  Matrix
 * @param  d  Dimension
 * @return    Sum of M along dimension d
 */
template <class T> inline static Matrix<T>
prod (const Matrix<T>& M, size_t d) {

	Matrix<size_t> sz = size(M);
	size_t        dim = sz[d];
	Matrix<T>     res;

	assert (d < M.NDim());

	// No meaningful sum over particular dimension
	if (dim == 1)
		return res;

	// Empty? allocation
	if (isempty(M))
		return res;

	// Inner and outer sizes
    Vector<size_t>::const_iterator ci = sz.Begin();
	size_t insize = std::accumulate (ci, ci+d, 1, c_multiply<size_t>);
    size_t outsize = std::accumulate (ci+d+1, ci+d+MIN(M.NDim(),sz.Size()), 1, c_multiply<size_t>);
        
    // Adjust size vector and allocate
	sz [d] = 1;
	res = Matrix<T>(sz);

	// Sum
#pragma omp parallel default (shared)
	{

#pragma omp for

		for (size_t i = 0; i < outsize; ++i) {

			for (size_t j = 0; j < insize; ++j) {
				res[i*insize + j] = T(0);
				for (size_t k = 0; k < dim; ++k)
					res[i*insize + j] += M[i*insize*dim + j + k*insize];
			}

		}

	}

	return res;

}


/**
 * @brief     Product of all elements
 *
 * Usage:
 * @code
 *   Matrix<cxfl> m   = rand<double> (8,7,6);
 *   m = prod (m);
 * @endcode
 *
 * @param  M  Matrix
 * @return    Sum of M along dimension d
 */
template <class T> inline static T
prod (const Matrix<T>& M) {
	return std::accumulate(M.Begin(), M.End(), T(1), c_multiply<T>);
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
SOS (const Matrix<T>& M, size_t d) {
	
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
	
	Vector<size_t> dim;
	Vector<float>  res;

	for (size_t i = 0; i < M.NDim(); ++i)
		if (size(M, i) > 1) {
			dim.push_back(size (M,i));
			res.push_back(resol(M,i));
		}
	
	Matrix<T> ret (dim);
	ret.Container() = M.Container();
	
	return ret;
	
}


/**
 * @brief           MATLAB-like permute
 *
 * Usage:
 * @code
 *   Matrix<cxfl> m   = rand<double> (2,3,4);
 *   m = permute (m, 0, 1, 2); // new dims: (4,2,3);
 * @endcode
 *
 * @param   M       Input matrix
 * @param   perm    New permuted dimensions
 * @return          Permuted matrix
 */

template<class T> inline static Matrix<T>
permute (const Matrix<T>& M, const size_t& n0, const size_t& n1, const size_t& n2) {
	Vector<size_t> odims = size(M);
	assert (numel(odims)==3); // Must be 3d
	Matrix<T> ret(odims[n0], odims[n1], odims[n2]);

	if        (n0 == 0) {
		if (n1 == 1) {
			assert (n2 == 2); // 0,1,2 : nothing to do
			return M;
		} else {
			assert (n2 == 1); // 0,2,1 : n1*n2 copies
			for (size_t k = 0; k < odims[n2]; ++k)
				for (size_t j = 0; j < odims[n1]; ++j)
					std::copy (&M(0,k,j), &M(0,k,j)+odims[0], &ret(0,j,k));
		}
	} else if (n0 == 1) {
		if (n1 == 0) {
			assert (n2 == 2); // 1,0,2
			for (size_t k = 0; k < odims[n2]; ++k)
				for (size_t j = 0; j < odims[n1]; ++j)
					for (size_t i = 0; i < odims[n0]; ++i)
						ret(i,j,k) = M(j,i,k);
		} else {
			assert (n2 == 0); // 1,2,0
			for (size_t k = 0; k < odims[n2]; ++k)
				for (size_t j = 0; j < odims[n1]; ++j)
					for (size_t i = 0; i < odims[n0]; ++i)
						ret(i,j,k) = M(k,i,j);
		}
	} else if (n0 == 2) {
		if (n1 == 1) {
			assert (n2 == 0); // 2,1,0
			for (size_t k = 0; k < odims[n2]; ++k)
				for (size_t j = 0; j < odims[n1]; ++j)
					for (size_t i = 0; i < odims[n0]; ++i)
						ret(i,j,k) = M(k,j,i);
		} else {
			assert (n2 == 1); // 2,0,1
			for (size_t k = 0; k < odims[n2]; ++k)
				for (size_t j = 0; j < odims[n1]; ++j)
					for (size_t i = 0; i < odims[n0]; ++i)
						ret(i,j,k) = M(j,k,i);
		}
	} else {

		assert(false);
	}

	return ret;

}



/**
 * @brief           MATLAB-like permute
 * 
 * Usage:
 * @code
 *   Matrix<cxfl> m   = rand<double> (1,8,7,1,6);
 *   m = permute (m); // dims: (8,7,6);
 * @endcode
 *
 * @param   M       Input matrix 
 * @param   perm    New permuted dimensions
 * @return          Permuted matrix
 */
template <class T> inline static Matrix<T>
permute (const Matrix<T>& M, const Matrix<size_t>& perm) {
	
	// Check that perm only includes one number between 0 and INVALID_DIM once
	size_t ndnew = perm.Size(), i, j;
	size_t ndold = ndims (M); 

	// Must have same number of dimensions
	assert (ndnew == ndold);

	// Every number between 0 and ndnew must appear exactly once
	Vector<bool> occupied;
	occupied.resize(ndnew);
	for (i = 0; i < ndnew; ++i) {
		assert (!occupied[perm[i]]);
		occupied [perm[i]] = true;
	}			

	// Old and new sizes
	Matrix<size_t> so = size (M);
	Matrix<size_t> sn (16,1);

	for (i = 0; i < ndnew; ++i)
		sn[i] = so[perm[i]];
	
	// Allocate new matrix with permuted dimensions
	Matrix<T> res(sn);

	// Relation of old to new indices
	size_t  d[16];
	size_t od[16];
	for (i = 0; i < ndnew; ++i) od[i] = perm[i];
	for (     ; i <    16; ++i)	od[i] =      i;
	
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
 * @brief          FLip up down
 * 
 * @param   M      Matrix
 * @return         Flipped matrix
 */

template <class T> inline static Matrix<T>
flipud (const Matrix<T>& M)  {

	size_t scol = size(M,0);
	size_t ncol = numel (M)/scol;

	Matrix<T> res = M;

    if (scol == 1) // trivial
        return res;

    typedef typename Vector<T>::iterator VI;
    
    VI rb = res.Container().begin();

	for (size_t i = 0; i < ncol; ++i)
        std::reverse(rb+i*scol, rb+(i+1)*scol);

	return res;

}

/**
 * @brief          FLip up down
 * 
 * @param  M        Matrix
 * @return         Flipped matrix
 */

template <class T> inline static Matrix<T>
fliplr (const Matrix<T>& M)  {

	size_t srow = size(M,1);
	size_t scol = size(M,0);
	size_t nrow = numel (M)/srow;

	Matrix<T> res (M.Dim());

    for (size_t i = 0; i < nrow; ++i)
        for (size_t j = 0; j < srow; ++j)
            res[j*scol+i] = M[(srow-1-j)*scol+i]; 

	return res;

}

#endif 
