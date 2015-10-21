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
	size_t rows = size(y,0);
    Matrix<T> ret(size(y));
    Matrix<RT> x = linspace<RT>(0,1,rows);
    Loess<T> l;
    Vector<T> yy;
    for (size_t i = 0; i < numel(y)/rows; ++i) {
    	yy = l.lowess(&x.Container()[0],&y.Container()[i*rows],rows);
        std::copy(yy.begin(), yy.end(), &ret.Container()[i*rows]);
    }
	return ret;
}

#endif /* SRC_MATRIX_INTERP_SMOOTH_HPP_ */
