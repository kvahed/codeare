/*
 * Smooth.hpp
 *
 *  Created on: Sep 29, 2015
 *      Author: kvahed
 */

#ifndef SRC_MATRIX_INTERP_SMOOTH_HPP_
#define SRC_MATRIX_INTERP_SMOOTH_HPP_

#include "TypeTraits.hpp"
#include "Matrix.hpp"
#include "PolyVal.hpp"
#include "loess.h"
#include <algorithm>
#include <cmath>
#include <fstream>

using namespace std;

template<class T> inline static Matrix<T>
smooth (const Matrix<T>& y, const size_t& span = 3, const INTERP::Method& method = INTERP::AKIMA) {
	typedef typename TypeTraits<T>::RT RT;
    Loess<T> l ((RT)span/(RT)size(y,0));
    Matrix<T> ret = l.lowess(y);
	return ret;
}

template<class T> inline static Matrix<T>
smooth (const View<T,false>& y, const size_t& span = 3, const INTERP::Method& method = INTERP::AKIMA) {
	Matrix<T> yy = y;
	return smooth(yy, span, method);
}

#endif /* SRC_MATRIX_INTERP_SMOOTH_HPP_ */
