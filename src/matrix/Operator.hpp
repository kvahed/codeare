/*
 * Operator.hpp
 *
 *  Created on: Apr 7, 2015
 *      Author: kvahed
 */

#ifndef SRC_MATRIX_OPERATOR_HPP_
#define SRC_MATRIX_OPERATOR_HPP_

#include "Matrix.hpp"

template<class T> class Operator {
public:
	Operator () {}
	virtual ~Operator () {}
	virtual Matrix<T> operator* (const Matrix<T>&) const {};
	virtual Matrix<T> operator->* (const Matrix<T>&) const {};
	virtual Matrix<T> operator/ (const Matrix<T>&) const {};
};


#endif /* SRC_MATRIX_OPERATOR_HPP_ */
