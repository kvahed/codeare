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
    typedef typename TypeTraits<T>::RT RT;
	Operator () {}
	virtual ~Operator () {}
	virtual Matrix<T> operator* (const MatrixType<T>&) const { return Matrix<T>(); }
	virtual Matrix<T> operator->* (const Matrix<T>&) const { return Matrix<T>(); }
	virtual Matrix<T> operator/ (const MatrixType<T>&) const { return Matrix<T>(); }
    virtual RT obj ( const Matrix<T>& x, const Matrix<T>& dx, const RT& t, RT& rmse) {return 0.;}
    virtual Matrix<T> df (const Matrix<T>& x) {return Matrix<T>();}
    virtual void Update (const Matrix<T>& dx) {}

};


#endif /* SRC_MATRIX_OPERATOR_HPP_ */
