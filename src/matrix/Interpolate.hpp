#ifndef __INTERPOLATE_HPP__
#define __INTERPOLATE_HPP__

#include "PolyVal.hpp"
#include "Access.hpp"
#include "IO.hpp" 

template <class T> inline static Matrix<T>
interp1 (Matrix<double>& x, Matrix<T>& y, const Matrix<double>& xi, const INTERP::Method& intm = INTERP::CSPLINE) {

	size_t  nxi = size(xi,0);
	size_t  nd  = size( y,1);
	size_t  nx  = size( x,0);

	Matrix<T> yi (nxi,nd);
	for (size_t j = 0; j < nd; j++) {
		
		PolyVal pv = PolyVal (x, (double*)&y[j*nx], intm);
		
		for (size_t i = 0; i < nxi; i++)
			yi [j * nxi + i] = pv.Lookup (xi[i]); 
		
	}
	
	return yi;
	
} 

#endif /* __INTERPOLATE_HPP__ */
