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

#include "Lapack.hpp"

template <class T> Matrix<T> 
Matrix<T>::prodt (const Matrix<T> &M) {
	
	return Lapack::GEMM (*this, M, 'C');
	
}


template <class T> Matrix<T> 
Matrix<T>::prod (const Matrix<T> &M, const char transa, const char transb) {
	
	return Lapack::GEMM (*this, M, transa, transb);
	
}


template<class T> T
Matrix<T>::Norm () const {
	
	T   res   = T(0.0);
	
	int n    = (int) Size();
	int incx = 1;
	
	if      (typeid(T) == typeid(  cxfl)) res = cblas_scnrm2 (n, &_M[0], incx);
	else if (typeid(T) == typeid(  cxdb)) res = cblas_zcnrm2 (n, &_M[0], incx);
	else if (typeid(T) == typeid(double)) res = cblas_dnrm2  (n, &_M[0], incx);
	else if (typeid(T) == typeid( float)) res = cblas_snrm2  (n, &_M[0], incx);
	
	else {
		for (int i = 0; i < Size(); i++)
			res += pow(_M[i],2);
		sqrt (res);
	}
	
	return res;
	
}


template<class T>  T 
Matrix<T>::dotc (Matrix<T>& M) const {
	
	T   res  = T(0.0);
	
	int n    = (int) Size();
	int incx = 1;
	
	if (typeid(T) == typeid(cxfl))
		cblas_cdotc_sub (n, &_M[0], incx, &M[0], incx, &res);
	else if (typeid(T) == typeid(cxdb))
		cblas_zdotc_sub (n, &_M[0], incx, &M[0], incx, &res);
	
	return res;
	
}

template<class T>  T 
Matrix<T>::dotu (Matrix<T>& M) const {
	
	T   res  = T(0.0);
	
	int n    = (int) Size();
	int incx = 1;
	
	if (typeid(T) == typeid(cxfl))
		cblas_cdotu_sub (n, &_M[0], incx, &M[0], incx, &res);
	else if (typeid(T) == typeid(cxdb))
		cblas_zdotu_sub (n, &_M[0], incx, &M[0], incx, &res);
	
	return res;
	
}

template<class T>  T 
Matrix<T>::dot (Matrix<T>& M) const {
	
	T   res  = T(0.0);
	
	int n    = (int) Size();
	int incx = 1;
	
	if (typeid(T) == typeid(cxfl))
	  cblas_cdotu_sub (n, &_M[0], incx, &M[0], incx, &res);
	else if (typeid(T) == typeid(cxdb))
	  cblas_zdotu_sub (n, &_M[0], incx, &M[0], incx, &res);
	else if (typeid(T) == typeid(double))
	  cblas_ddot_sub (n, &_M[0], incx, &M[0], incx, &res);
	else if (typeid(T) == typeid(float))
	  cblas_sdot_sub (n, &_M[0], incx, &M[0], incx, &res);
	
	return res;
	
}

