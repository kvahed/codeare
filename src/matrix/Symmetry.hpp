/*
 * Symmtry.hpp
 *
 *  Created on: Oct 23, 2015
 *      Author: kvahed
 */

#ifndef SRC_MATRIX_SYMMTRY_HPP_
#define SRC_MATRIX_SYMMTRY_HPP_

#include "Algos.hpp"

template<class T> inline static bool issymmetric (const Matrix<T>& A) {
	if (!is2d(A))
		return false;
	if (TypeTraits<T>::IsComplex())
		return false;
	size_t m = size(A,0);
	for (size_t j = 0; j < m; ++j)
		for (size_t i = 0; i < m; ++i)
			if (A(i,j)!=A(j,i))
				return false;
	return true;
}
template<class T> inline static bool ishermitian (const Matrix<T>& A) {
	if (!is2d(A))
		return false;
	if (TypeTraits<T>::IsReal())
		return false;
	size_t m = size(A,0);
	for (size_t j = 0; j < m; ++j)
		for (size_t i = 0; i < m; ++i)
			if (A(i,j)!=TypeTraits<T>::Conj(A(j,i)))
				return false;
	return true;
}


#endif /* SRC_MATRIX_SYMMTRY_HPP_ */
