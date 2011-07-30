#include <cstdlib>
#include <cblas.h>

template <class T>
Matrix<T> 
Matrix<T>::prodt (Matrix<T> &M) {
	
	M = M.tr();
	return this->prod(M);

}


template <class T>
Matrix<T> 
Matrix<T>::prod (Matrix<T> &M) {
	
    assert(Dim(1) == M.Dim(0));
	return GEMM(M, 'N');
	
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

	Matrix<T> res (m,n);

	if (typeid(T) == typeid(double))
		dgemm_ (&transa, &transb, &m, &n, &k, &alpha, &at(0), &lda, &M.at(0), &ldb, &beta, &res.at(0), &ldc);
	else if (typeid(T) == typeid(raw))
		cgemm_ (&transa, &transb, &m, &n, &k, &alpha, &at(0), &lda, &M.at(0), &ldb, &beta, &res.at(0), &ldc);

	return res;
	
}



template<class T>
T
Matrix<T>::norm () const {

	T   res   = (T)0.0;

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
	
	T   res  = (T)0.0;

	int n    = Size();
	int incx = 1;

	if (typeid(T) == typeid(raw))
		cblas_cdotc_sub (n, &_M[0], incx, &M[0], incx, &res);
	
	return res;

}
