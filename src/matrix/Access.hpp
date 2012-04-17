#ifndef __ACCESS_HPP__
#define __ACCESS_HPP__

#include "Algos.hpp"

/**
 * @brief              Get a volume
 * 
 * @param  M           Matrix
 * @param  s           # of volume
 * @return             Desired volume of M
 */
template <class T> static inline Matrix<T> 
Volume (const Matrix<T>& M, const size_t s) {
	
	Matrix<T> res (size(M,0), size(M,1), size(M,2));
	size_t nc = numel (res);
	
	memcpy (&res[0], M.Data(s*nc), nc * sizeof(T));
	
	return res;
	
}
	

/**
 * @brief              Set a volume in a matrix
 * 
 * @param  M           Matrix to insert to
 * @param  s           # of volume
 * @param  A           Matrix to insert
 */
template <class T> inline static void
Volume (Matrix<T>& M, const size_t s, const Matrix<T> A) {
	
	assert (size(M,0) == size(A,0));
	assert (size(M,1) == size(A,1));
	assert (size(M,2) == size(A,2));

	size_t nc = size(M,0) * size(M,1) * size(M,2);

	memcpy (&M[s*nc], A.Data(), nc*sizeof(T));
	
}
	

/**
 * @brief              Get a slice
 * 
 * @param  M           Matrix
 * @param  s           # of slice
 * @return             Desired slice of M
 */
template <class T> static inline Matrix<T> 
Slice (const Matrix<T>& M, const size_t s) {
	
	Matrix<T> res (size(M,0),size(M,1));
	size_t nc = numel (res);
	
	memcpy (&res[0], M.Data(s*nc), nc*sizeof(T));
	
	return res;
	
}

	
/**
 * @brief              Set a slice
 * 
 * @param  M           Matrix
 * @param  s           # of slice
 * @param  A           New slice
 */
template <class T> static inline void
Slice (Matrix<T>& M, const size_t& s, const Matrix<T> A) {
	
	assert (size(M,0) == size(A,0));
	assert (size(M,1) == size(A,1));

	size_t ns = size(M,0) * size(M,1);

	memcpy (&M[s * ns], A.Data(), ns * sizeof(T));
	
}

	
/**
 * @brief              Set a slice
 * 
 * @param  M           Matrix
 * @param  s           # of slice
 * @param  v           Scalar value           
 */
template <class T> static inline void
Slice (Matrix<T>& M, const size_t& s, const T& v) {
	
	size_t ns = size(M,0) * size(M,1);

	vector<T> vv; 
	vv.resize(ns, v);

	memcpy (&M[s * ns], &vv[0], ns * sizeof(T));
	
}

	
/**
 * @brief              Get a row
 * 
 * @param  M           Matrix
 * @param  r           # of row
 * @return             Desired row of M
 */
template <class T> static inline Matrix<T> 
Row (const Matrix<T>& M, const size_t r)  {
	
	Matrix<T> res (size(M, 1),1);
	
	for (size_t i = 0; i < size(M, 1); i++)
		res[i] = M[r + i*size(M, 0)];
	
	return res;
	
}


/**
 * @brief              Copy a row (A) into matrix M
 * 
 * @param  M           Matrix to copy into
 * @param  r           # of row
 * @param  A           Matrix to copy from
 */
template <class T> static inline void 
Row (Matrix<T>& M, const size_t r, const Matrix<T> A)  {
	
	size_t nc, nr, ns, s, rr;

	nc = size(M,1);
	nr = size(M,0);
	ns = nc * nr;

	assert (size(M,1) == numel (A));
	
	s  = floor (r/nr);
	rr = r % nr;

	for (size_t i = 0; i < nc; i++)
		M[rr + i * nr + s * ns] = A[i];
	
}
	
	
/**
 * @brief              Get a row
 * 
 * @param  M           Matrix
 * @param  c           # of row
 * @return             Desired row of M
 */
template <class T> static inline Matrix<T> 
Column (const Matrix<T>& M, const size_t c) {
	
	Matrix<T> res (size(M, 0),1);
	
	memcpy (&res[0], M.Data(c*size(M, 0)), size(M, 0) * sizeof(T));
	
	return res;
	
}


/**
 * @brief              Get a row
 * 
 * @param  M           Matrix
 * @param  c           # of row
 * @param  A           Vector to copy from
 */
template <class T> static inline void
Column (Matrix<T>& M, const size_t c, const Matrix<T> A) {
	
	assert (size(M,0) == size (A,0));
	size_t nc = size(M,0);

	memcpy (&M[c*nc], A.Data(), nc * sizeof(T));
	
}


#endif
