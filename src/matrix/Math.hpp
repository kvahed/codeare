#include "Matrix.hpp"


/**
 * @brief          Exponential
 * 
 * @param  m       Matrix
 * @return         Exponential
 */ 
template <class T> inline Matrix<T>
exp (const Matrix<T>& m) {
	
	Matrix<T> res = m;
	
	for (size_t i = 0; i < res.Size(); i++)
		res [i] = exp (res[i]);
	
	return res;
	
}



/**
 * @brief          Logarithm
 * 
 * @param  m       Matrix
 * @return         Logarithm
 */ 
template <class T> inline Matrix<T>
log (const Matrix<T>& m) {
	
	Matrix<T> res = m;
	
	for (size_t i = 0; i < res.Size(); i++)
		res [i] = (res[i] != 0) ? log (res[i]) : 0.0;
	
	return res;
	
}



/**
 * @brief          Sine
 * 
 * @param  m       Matrix
 * @return         Sine
 */ 
template <class T> inline Matrix<T>
sin (const Matrix<T>& m) {
	
	Matrix<T> res = m;
	
	for (size_t i = 0; i < res.Size(); i++)
		res [i] = sin (res[i]);
	
	return res;
	
}



/**
 * @brief          Cosine
 * 
 * @param  m       Matrix
 * @return         Cosine
 */ 
template <class T> inline Matrix<T>
cos (const Matrix<T>& m) {
	
	Matrix<T> res = m;
	
	for (size_t i = 0; i < res.Size(); i++)
		res [i] = cos (res[i]);
	
	return res;
	
}


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
