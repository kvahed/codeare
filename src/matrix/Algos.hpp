
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

#define NOMINMAX

#include <Matrix.hpp>
#include <math.h>
#include <limits>
#include <vector>
#include <algorithm>    // std::reverse
#include <numeric>


#if !defined(_MSC_VER) || _MSC_VER>1200
#include <boost/math/special_functions/fpclassify.hpp>

template<class T> inline static bool is_nan (T const& x) {
    return boost::math::isnan (x);
}
template <class T> inline static bool is_inf (T const& x) {
    return boost::math::isinf (x);
}
#endif


template<class T, class S>
inline static bool eq (const MatrixType<T>& A, const MatrixType<S>& B) {
    assert (A.Size() == B.Size());
    for (size_t i = 0; i < A.Size(); ++i)
        if (A[i]!=B[i])
            return false;
    return true;
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
template <class T> inline static  size_t nnz (const Matrix<T>& M) {
	size_t nz   = 0;
	for (size_t i = 0; i < M.Size(); ++i)
		if (M[i] != T(0))
			++nz;
	return nz;
}

template<class T, class S> inline unsigned short issame (const Matrix<T>& A, const Matrix<S>& B) {
	if (numel(A) != numel(B))
		return 0;
	for (auto i = 0; i < numel(A); ++i)
		if (A[i]!=B[i])
			return 0;
	if (size(A)==size(B))
		return 2;
	return 1;
}


/**
 * @brief     Is matrix X-dimensional?
 *
 * @param  M  Matrix
 * @param  d  Dimension
 * @return    X-dimensional?
 */
template <class T>  inline static  bool isxd (const Matrix<T>& M, size_t d) {

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
template <class T>  inline static bool isvec (const Matrix<T>& M) {
	
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
isempty (const MatrixType<T>& M) {
	
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
template <class T> inline static size_t ndims (const MatrixType<T>& M) {
	
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
template <class T> inline static  Matrix<T> diag (const Matrix<T>& M) {
	assert (is2d(M));
	size_t sz = (std::min)(size(M,0),size(M,1));
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
template <class T, paradigm P> inline static size_t numel (const MatrixType<T,P>& M) {
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
size               (const MatrixType<T>& M, size_t d) {
	return M.Dim(d);
}
template <class T>  size_t
size               (const MatrixType<T,MPI>& M, size_t d) {
	return M.Dim(d);
}




/**
 * @brief           Get vector of dimensions
 *
 * @param   M       Matrix
 * @return          Dimension vector.
 */
template <class T,paradigm P>  inline static  Vector<size_t>
size               (const MatrixType<T,P>& M) {
	return M.Dim();
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
#ifdef _MSC_VER
#  ifdef max
#    undef max
#  endif
#endif
template<class T> inline static Matrix<T> max (const Matrix<T>& M, const size_t& dim = 0) {
	Vector<size_t> dims = size(M); size_t m = dims[0]; size_t n = numel(M)/m;
	dims[dim] = 1;
	Matrix<T> ret(dims);
	for (size_t i = 0; i < n; ++i)
		ret[i] = *std::max_element(M.Begin()+i*m,M.Begin()+(i+1)*m);
	return ret;
}
template<class T> inline static Matrix<T> max (const View<T,true>& M, const size_t& dim = 0) {
	Vector<size_t> dims = size(M); size_t m = dims[0]; size_t n = numel(M)/m;
	dims.erase(dims.begin());
	Matrix<T> ret(dims);
	for (size_t j = 0; j < n; ++j) {
		ret[j] = -1e20;
		for (size_t i = 0; i < m; ++i)
			if (ret[j]<=M[j*m+i]) ret[j] = M[j*m+i];
	}
	return ret;
}
template<class T> inline static T mmax (const Matrix<T>& M) {
	return *std::max_element(M.Begin(), M.End());
}
template <class T> inline static T mmax (const View<T, true>& V) {
	T mx = -1e20;
	for (size_t i = 0; i < numel(V); ++i)
		if (V[i] >= mx)
			mx = V[i];
	return mx;
}

/**
 * @brief           Maximal element
 *
 * @param  M        Matrix
 * @return          Maximum
 */
#ifdef _MSC_VER
#  ifdef min
#    undef min
#  endif
#endif
template<class T> inline static Matrix<T> min (const Matrix<T>& M, const size_t& dim = 0) {
	Vector<size_t> dims = size(M); size_t m = dims[0]; size_t n = numel(M)/m;
	dims[dim] = 1;
	Matrix<T> ret(dims);
	for (size_t i = 0; i < n; ++i)
		ret[i] = *std::min_element(M.Begin()+i*m,M.Begin()+(i+1)*m);
	return ret;
}
template<class T> inline static Matrix<T> min (const View<T,true>& M, const size_t& dim = 0) {
	Vector<size_t> dims = size(M); size_t m = dims[0]; size_t n = numel(M)/m;
	dims.erase(dims.begin());
	Matrix<T> ret(dims);
	for (size_t j = 0; j < n; ++j) {
		ret[j] = 1e20;
		for (size_t i = 0; i < m; ++i)
			if (ret[j]>=M[j*m+i]) ret[j] = M[j*m+i];
	}
	return ret;
}
template<class T> inline static T mmin (const Matrix<T>& M) {
	return *std::min_element(M.Begin(), M.End());
}
template <class T> inline static T mmin (const View<T, true>& V) {
	T mx = 1e-20;
	for (size_t i = 0; i < numel(V); ++i)
		if (V[i] <= mx)
			mx = V[i];
	return mx;
}

#include "CX.hpp"
/**
 * @brief           Transpose
 *
 * @param  M        2D Matrix
 * @param  c        Conjugate while transposing
 *
 * @return          Non conjugate transpose
 */
template <class T> inline static  Matrix<T> transpose (const Matrix<T>& M, bool c = false) {
	assert (is2d(M));
	Matrix<T> res (size(M,1),size(M,0));
	for (size_t j = 0; j < size(res,1); ++j)
		for (size_t i = 0; i < size(res,0); ++i)
			res(i,j) = M(j,i);
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
template <class T> inline static  Matrix<T> resize (const Matrix<T>& M, size_t sz) {

	Matrix<T> res (sz,1);
	size_t copysz = std::min(numel(M), sz);

    typename Vector<T>::      iterator rb = res.Begin ();
    typename Vector<T>::const_iterator mb =   M.Begin ();
    
    std::copy (mb, mb+copysz, rb);

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
template <class T> inline static  Matrix<T> resize (const Matrix<T>& M, size_t sc, size_t sl) {
	assert(sl*sc==numel(M));
	Matrix<T> ret(sc,sl);
	ret.Container() = M.Container();
	return ret;
}


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
template <class T> inline static Matrix<T> resize (const Matrix<T>& M, const size_t& s0,
		const size_t& s1, const size_t& s2) {
	assert (numel(M)==s0*s1*s2);
	Matrix<T> res (s0,s1,s2);
	res.Container() = M.Container();
	return res;
}

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
template <class T> inline static Matrix<T> resize (const Matrix<T>& M, const size_t& s0,
		const size_t& s1, const size_t& s2, const size_t& s3) {
	assert (numel(M)==s0*s1*s2*s3);
	Matrix<T> res (s0,s1,s2,s3);
	res.Container() = M.Container();
	return res;
}

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
template <class T> inline static Matrix<T> resize (const Matrix<T>& M, const size_t& s0,
		const size_t& s1, const size_t& s2, const size_t& s3, const size_t& s4) {
	assert (numel(M)==s0*s1*s2*s3*s4);
	Matrix<T> res (s0,s1,s2,s3,s4);
	res.Container() = M.Container();
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
	size_t copysz  = std::min(numel(M), numel(res));

    typename Vector<T>::      iterator rb = res.Begin ();
    typename Vector<T>::const_iterator mb =   M.Begin ();

    std::copy (mb, mb+copysz, rb);

	return res;
	
}

template <class T> inline static T sum2 (const Matrix<T>& M) {
	return std::accumulate (M.Begin(), M.End(), (T)1, std::plus<T>());
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
template <class T> inline static Matrix<T> sum (const MatrixType<T>& M, const size_t& d = 0) {
	
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
	for (size_t i = d+1; i < std::min(M.NDim(),sz.size()); ++i)
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


template <class T> inline static Matrix<T> mean (const MatrixType<T>& M, const size_t& d = 0) {
	return sum(M,d)/size(M,d);
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
template <class T> inline static Matrix<T> prod (const Matrix<T>& M, size_t d) {

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
    size_t outsize = std::accumulate (ci+d+1, ci+d+std::min(M.NDim(),sz.Size()), 1, c_multiply<size_t>);
        
    // Adjust size vector and allocate
	sz [d] = 1;
	res = Matrix<T>(sz);

	// Sum
#pragma omp parallel default (shared)
	{
#pragma omp for
		for (size_t i = 0; i < outsize; ++i)
			for (size_t j = 0; j < insize; ++j) {
				res[i*insize + j] = T(0);
				for (size_t k = 0; k < dim; ++k)
					res[i*insize + j] += M[i*insize*dim + j + k*insize];
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
template <class T> inline static T prod (const Matrix<T>& M) {
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
template <class T> inline static Matrix<typename TypeTraits<T>::RT>
    sos (const Matrix<T>& M, long d = -1) {
    typedef typename TypeTraits<T>::RT real_type;
    if (d == -1)
        d = ndims(M)-1;
    assert (d <= ndims(M)-1);
    const real_type* rt = (real_type*)&M[0];
	Matrix<real_type> res(size(M));
    for (size_t i = 0; i < numel(M); ++i)
        res[i] = rt[2*i]*rt[2*i]+rt[2*i+1]*rt[2*i+1];
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
template <class T> inline static Matrix<T> squeeze (const Matrix<T>& M) {
    Matrix<T> ret = M;
    ret.Squeeze();
	return ret;
}
template<class T> inline static Matrix<T> squeeze (const View<T,true>& V) {
	Vector<size_t> vdim = size(V), dim;
	for (size_t i = 0; i < vdim.size(); ++i)
		if (vdim[i] > 1)
			dim.push_back(vdim[i]);
    Matrix<T> ret(dim);
    for (size_t i = 0; i < numel(V); ++i)
        ret[i] = V[i];
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

template<class T> inline static Matrix<T> permute (const Matrix<T>& M, const size_t& n0,
		const size_t& n1, const size_t& n2) {
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
 *   Matrix<cxfl> m   = rand<double> (2,3,4);
 *   m = permute (m, 0, 1, 2); // new dims: (4,2,3);
 * @endcode
 *
 * @param   M       Input matrix
 * @param   perm    New permuted dimensions
 * @return          Permuted matrix
 */

template<class T> inline static Matrix<T> permute (const Matrix<T>& M, const size_t& n0,
		const size_t& n1) {
	Vector<size_t> odims = size(M);
	assert (numel(odims)==2); // Must be 3d
	Matrix<T> ret(odims[n0], odims[n1]);

	if        (n0 == 0) {// 0,1: nothing to do
		assert (n1 == 1);
		return M;
	} else if (n0 == 1) {// 1,0: transpose
		assert (n1 == 0);
		for (size_t j = 0; j < odims[n1]; ++j)
			for (size_t i = 0; i < odims[n0]; ++i)
				ret(i,j) = M(j,i);
	} else { // We should never be here
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
//#include "Print.hpp"
template <class T> inline static Matrix<T> permute (const Matrix<T>& M, const Vector<size_t>& perm) {
	
	// Check that perm only includes one number between 0 and INVALID_DIM once
	size_t ndnew = perm.size(), i = 0;
	size_t ndold = ndims (M); 

	// Must have same number of dimensions
	assert (ndnew == ndold);

	// Every number between 0 and ndnew must appear exactly once
	Vector<cbool> occupied;
	occupied.resize(ndnew);
	for (i = 0; i < ndnew; ++i) {
		assert (!occupied[perm[i]]);
		occupied [perm[i]] = true;
	}			

	// Old and new sizes
	Vector<size_t> so = size (M);
	Vector<size_t> sn (16,1);

	for (i = 0; i < ndnew; ++i)
		sn[i] = so[perm[i]];
	
	// Allocate new matrix with permuted dimensions
	Matrix<T> res(sn);

	// Relation of old to new indices
    Vector<size_t>  d(16), od(16);
	for (i = 0; i < ndnew; ++i) od[i] = perm[i];
	for (     ; i <    16; ++i)	od[i] =      i;

	// Copy data accordingly
	for (d[15] = 0; d[15] < size(M,15); d[15]++)
		for (d[14] = 0; d[14] < size(M,14); d[14]++)
			for (d[13] = 0; d[13] < size(M,13); d[13]++)
				for (d[12] = 0; d[12] < size(M,12); d[12]++)
					for (d[11] = 0; d[11] < size(M,11); d[11]++)
						for (d[10] = 0; d[10] < size(M,10); d[10]++)
							for (d[ 9] = 0; d[ 9] < size(M, 9); d[ 9]++)
								for (d[ 8] = 0; d[ 8] < size(M, 8); d[ 8]++)
									for (d[ 7] = 0; d[ 7] < size(M, 7); d[ 7]++)
										for (d[ 6] = 0; d[ 6] < size(M, 6); d[ 6]++)
											for (d[ 5] = 0; d[ 5] < size(M, 5); d[ 5]++)
												for (d[ 4] = 0; d[ 4] < size(M, 4); d[ 4]++)
													for (d[ 3] = 0; d[ 3] < size(M, 3); d[ 3]++)
														for (d[ 2] = 0; d[ 2] < size(M, 2); d[ 2]++)
															for (d[ 1] = 0; d[ 1] < size(M, 1); d[ 1]++)
																for (d[ 0] = 0; d[ 0] < size(M, 0); d[ 0]++) 
																	res (d[od[ 0]],d[od[ 1]],d[od[ 2]],d[od[ 3]], 
																		 d[od[ 4]],d[od[ 5]],d[od[ 6]],d[od[ 7]], 
																		 d[od[ 8]],d[od[ 9]],d[od[10]],d[od[11]],
																		 d[od[12]],d[od[13]],d[od[14]],d[od[15]]) =
																	  M (d[ 0],d[ 1],d[ 2],d[ 3],d[ 4],d[ 5],
																		 d[ 6],d[ 7],d[ 8],d[ 9],d[10],d[11],
																		 d[12],d[13],d[14],d[15]);

    // Remove trailing singelton dimensions
    sn.erase(sn.begin()+ndnew, sn.begin()+16);
	return resize(res, sn);

																		
}


/**
 * @brief           MATLAB-like diff
 */
template<class T> inline static Matrix<T> diff (const Matrix<T>& rhs, const size_t& n = 1,
		const size_t& dim = 0) {
    Matrix<T> ret, tmp;
    Vector<size_t> rhssize = size(rhs), pdims = rhssize, dim_ord;
    size_t ndims = rhssize.size();

    if (dim > ndims) {
        printf ("  *** ERROR (%s:%d) - diff dimension %zu"
                "exceeds matrix dimension %zu", __FILE__, __LINE__, dim, ndims);
        throw 404;
    }

    if (dim == 0) {
        ret = rhs;
    } else { 
        if (dim == 1) {
            if (ndims == 2) {
                ret = permute (rhs,1,0);
            } else if (ndims == 3) {
                ret = permute (rhs,1,0,2);
            } else {
                dim_ord.resize(ndims);
                std::iota(dim_ord.begin(), dim_ord.end(), 1);
                std::iter_swap(dim_ord.begin(), dim_ord.begin()+1);
                ret = permute (rhs,dim_ord);
            }
        } else if (dim == 2) {
            if (ndims == 3)
                ret = permute (rhs,2,0,1);
            else {
                dim_ord.resize(ndims);
                std::iota(dim_ord.begin(), dim_ord.end(), 1);                
                std::iter_swap(pdims.begin(), pdims.begin()+2);
                ret = permute (rhs,dim_ord);
            }
        }
    }
    
    
    Vector<T> col (size(rhs,0));
    typename Vector<T>::iterator b, e, m;
    typename Vector<T>::const_iterator rb;
    size_t col_len = col.size();
    for (size_t i = 0 ; i < n; ++i) {
        tmp = ret;
        for (b = ret.Begin(), e = b+col_len, m = b+1, rb = tmp.Begin();
             b < ret.End(); b += col_len, e += col_len, m += col_len, rb += col_len) {
            std::rotate (b, m, e);
            std::transform (b, e, rb, b, std::minus<T>());
        }
    }

    
    if (dim == 1) {
        if (ndims == 2) {
            ret = permute (rhs,1,0);
        } else if (ndims == 3) {
            ret = permute (rhs,1,0,2);
        } else {
            dim_ord.resize(ndims);
            std::iota(dim_ord.begin(), dim_ord.end(), 1);
            std::iter_swap(dim_ord.begin(), dim_ord.begin()+1);
            ret = permute (rhs,dim_ord);
        }
    } else if (dim == 2) {
        if (ndims == 3)
                ret = permute (rhs,2,0,1);
        else {
            dim_ord.resize(ndims);
            std::iota(dim_ord.begin(), dim_ord.end(), 1);                
            std::iter_swap(pdims.begin(), pdims.begin()+2);
            ret = permute (rhs,dim_ord);
        }
    }
    
    return ret;
    
}



/**
 * @brief          FLip up down
 * 
 * @param   M      Matrix
 * @return         Flipped matrix
 */
template <class T> inline static Matrix<T> flipud (const Matrix<T>& M)  {

	size_t scol = size(M,0), ncol = numel(M)/scol;
	Matrix<T> res = M;

    if (scol == 1) // trivial
        return res;

    typedef typename Vector<T>::iterator VI;
	for (VI i = res.Container().begin(); i < res.Container().end(); i += scol)
        std::reverse(i, i+scol);
	return res;

}

/**
 * @brief          FLip left right
 * 
 * @param  M        Matrix
 * @return         Flipped matrix
 */
template <class T> inline static Matrix<T> fliplr (const Matrix<T>& M)  {

	size_t srow = size(M,1), scol = size(M,0), nrow = numel (M)/srow;
	Matrix<T> res (M.Dim());

    for (size_t i = 0; i < nrow; ++i)
        for (size_t j = 0; j < srow; ++j)
            res[j*scol+i] = M[(srow-1-j)*scol+i]; 

	return res;

}

/**
 * @brief Sort keep original indices
 */
typedef enum sort_dir {
	ASCENDING, DESCENDING
} sort_dir;

/**
 * @brief   Get sort indices sorting elements of m
 * @param  m Data to sort
 * @param  sd Sort direction
 * @return sort indices
 */
template <typename T> inline static Vector<size_t> sort (const Matrix<T> &m,
		const sort_dir sd = ASCENDING) {
	Vector<size_t> idx(m.Size());
	std::iota(idx.begin(), idx.end(), 0);
	if (sd == ASCENDING)
		sort(idx.begin(), idx.end(), [&m](size_t i1, size_t i2) {return m[i1] < m[i2];});
	else
		sort(idx.begin(), idx.end(), [&m](size_t i1, size_t i2) {return m[i1] > m[i2];});

	return idx;
}

#endif 
