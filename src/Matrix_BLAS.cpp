/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but 
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301  USA
 */

#include "BLAS.hpp"
#include <cstdlib>
#include <cblas.h>

template <class T>
Matrix<T> 
Matrix<T>::prodt (Matrix<T> &M) {
	
	Matrix<T> tmp = M.tr();
	return this->prod(tmp);
	
}


template <class T>
Matrix<T> 
Matrix<T>::prod (Matrix<T> &M) {
	
    assert (Dim(1) == M.Dim(0));
	return GEMM (M, 'N');
	
}


template<class T>
Matrix<T> 
Matrix<T>::GEMM (Matrix<T>& M, char transb) {
	
	char transa = 'N';
	
    int  m      =   Dim(0);
    int  n      = M.Dim(1);
    int  k      =   Dim(1);
	int  lda    =   m;
	int  ldb    =   k;
	int  ldc    =   m;
	
	T    alpha  =   T(1.0);
	T    beta   =   T(0.0);
	
	Matrix<T> res (m, (transb == 'N') ? M.Dim(1) : M.Dim(0));
	
	if (typeid(T) == typeid(double))
		dgemm_ (&transa, &transb, &m, &n, &k, &alpha, &_M[0], &lda, &M[0], &ldb, &beta, &res[0], &ldc);
	else if (typeid(T) == typeid(raw))
		cgemm_ (&transa, &transb, &m, &n, &k, &alpha, &_M[0], &lda, &M[0], &ldb, &beta, &res[0], &ldc);
	
	return res;
	
}



template<class T>
T
Matrix<T>::Norm () const {
	
	T   res   = T(0.0);
	
	int n    = Size();
	int incx = 1;
	
	if      (typeid(T) == typeid(   raw)) res = cblas_scnrm2 (n, _M, incx);
	else if (typeid(T) == typeid(double)) res = cblas_dnrm2  (n, _M, incx);
	
	else {
		for (int i = 0; i < Size(); i++)
			res += pow(_M[i],2);
		sqrt (res);
	}
	
	return res;
	
}


template<class T>
T 
Matrix<T>::dotc (Matrix<T>& M) const {
	
	T   res  = T(0.0);
	
	int n    = Size();
	int incx = 1;
	
	if (typeid(T) == typeid(raw))
		cblas_cdotc_sub (n, &_M[0], incx, &M[0], incx, &res);
	
	return res;
	
}
