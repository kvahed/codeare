#include "Matrix.hpp"

/**
 * @brief          Cumulative sum of all elements
 * 
 * @param  m       Vector
 * @return         Vector of cumulative sums
 */ 
template <class T> inline Matrix<T> 
cumsum (const Matrix<T>& m) {

	Matrix<T> res = m;

	for (size_t i = 1; i < res.Size(); i++)
		res [i] += res[i-1];

	return res;

}
