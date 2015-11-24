/*
 * LocalMaxima.hpp
 *
 *  Created on: Oct 20, 2015
 *      Author: kvahed
 */

#ifndef SRC_MATRIX_CURVES_LOCALMAXIMA_HPP_
#define SRC_MATRIX_CURVES_LOCALMAXIMA_HPP_


#include <Algos.hpp>
#include <Math.hpp>
#include <Print.hpp>
#include <climits>

template<class T> inline static Vector<size_t>
findLocalMaxima (const Matrix<T>& A, const size_t min_dist = 1, const T& min_val = -T(INT_MIN)) {
	Vector<size_t> ret;
	size_t m = size(A,0);
	Matrix<size_t> idx = linspace<size_t>(0,size(A,0)-1,size(A,0));
	Matrix<short> s = sign(diff(A));
	for (size_t i = 1; i < m; ++i)
		if (s[i-1]==1 && s[i]==-1 && A[i]>min_val) {
			ret.push_back(i);
            i += min_dist-1;
        }
	return ret;
}
#endif /* SRC_MATRIX_CURVES_LOCALMAXIMA_HPP_ */
