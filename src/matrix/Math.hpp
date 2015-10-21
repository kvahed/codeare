#include "Matrix.hpp"

/**
 * @brief          Cumulative sum of all elements
 * 
 * @param  m       Vector
 * @return         Vector of cumulative sums
 */ 
template <class T> inline Matrix<T> cumsum (const Matrix<T>& m) {
	Matrix<T> res = m;
	for (size_t i = 1; i < res.Size(); i++)
		res [i] += res[i-1];
	return res;
}

template<class T> inline static Matrix<T> exp (const Matrix<T>& M) {
    Matrix<T> ret(M.Dim());
    for (size_t i = 0; i < M.Size(); ++i)
        ret[i] = exp(M[i]);
    return M;
}

template<class T> inline Matrix<short> sign (const Matrix<T>& A) {
	Matrix<short> ret(size(A));
	for (size_t i=0; i < A.Size(); ++i)
		ret[i] = (A[i]>=0) ? 1 : -1;
	return ret;
}
