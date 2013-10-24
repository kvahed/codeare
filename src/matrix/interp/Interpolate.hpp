#ifndef __INTERPOLATE_HPP__
#define __INTERPOLATE_HPP__

#include "PolyVal.hpp"
#include "Algos.hpp"
#include "Access.hpp"

/**
 * @brief            1D Interpolation 
 *
 * @param  x         Original base
 * @param  y         Original values
 * @param  xi        Interpolation base
 * @param  intm      Interpolation method @see INTERP::Method
 * @return           Interpolation values
 */
template <class T> inline static Matrix<T>
interp1 (const Matrix<double>& x, const Matrix<T>& y,
		const Matrix<double>& xi, const INTERP::Method& intm = INTERP::CSPLINE) {

	size_t  nx = size(x,0);
	assert (nx > 0);
	assert (nx == size(y,0));

	size_t  nd  = numel(y)/nx;
	size_t  nxi  = size(xi,0);

	Matrix<T> yi (nxi,nd);
	for (size_t j = 0; j < nd; j++) {

		PolyVal<T> pv (x, (T*) y.Ptr(j*nx), intm);
		
		for (size_t i = 0; i < nxi; i++)
			yi [j * nxi + i] = pv.Lookup (xi[i]); 
		
	}
	
	return yi;
	
} 

#endif /* __INTERPOLATE_HPP__ */
