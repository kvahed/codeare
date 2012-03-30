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
	size_t nc = M.Dim(0)*M.Dim(1)*M.Dim(2);
	
	memcpy (&res[0], M.Data(s*nc), nc * sizeof(T));
	
	return res;
	
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
	size_t nc = M.Dim(0)*M.Dim(1);
	
	memcpy (&res[0], M.Data(s*nc), nc*sizeof(T));
	
	return res;
	
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
	
	Matrix<T> res (M.Dim(1),1);
	
	for (size_t i = 0; i < M.Dim(1); i++)
		res[i] = M[r + i*M.Dim(0)];
	
	return res;
	
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
	
	Matrix<T> res (M.Dim(0),1);
	
	memcpy (&res[0], &M[c*M.Dim(0)], M.Dim(0) * sizeof(T));
	
	return res;
	
}


#endif
