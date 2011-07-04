#include <cstdlib>
#include <cblas.h>

template <class T>
Matrix<T> 
Matrix<T>::prod(Matrix<T> &M) {

    assert(width() == M.height());

	if (typeid(T) == typeid(double) || typeid(T) == typeid(raw)) // Fast BLAS code
		return GEMM(M);

	else {                                                          // Standard impl

		Matrix<T> res;
		
		res.Dim(0)  = M.Dim(COL);//_dim[0];
		res.Dim(1)  = _dim[LIN];//M.Dim(1);
		res.Reset();
		
		for (int i = 0; i < res.height(); i++)
			for (int j = 0; j < res.width(); j++)
				for (int k = 0; k < width(); k++) {
					res[i * res.width() + j] += _M[i * width() + k] * M[k * M.width() + j];
					if (res[i * res.width() + j] < 0 || res[i * res.width() + j] > ICE_SHRT_MAX) {
						res[i * res.width() + j] = ICE_SHRT_MAX;
						break;
					}
				}


		
		return res;
		
	}
	
}

template<class T>
Matrix<T> 
Matrix<T>::GEMM (Matrix<T>& M) {
	
	Matrix<T> res;

	char transa = 'N';
	char transb = 'N';

	int i = 0, j = 0, l = 0;

	int  m      = _dim[LIN];
	int  n      = M.Dim(COL);

	int  k      = _dim[COL];

	T    alpha  = T(1);

	T*   a      = new T[Size()];
	i = 0;
	for (j = 0; j < _dim[COL]; j++)
		for (l = 0; l < _dim[LIN]; l++, i++)
			a[i] = _M[j+l*_dim[COL]];

	int  lda    = _dim[LIN];

	T*   b      = new T[M.Size()];
	i = 0;
	for (j = 0; j < M.Dim(COL); j++)
		for (l = 0; l < M.Dim(LIN); l++, i++)
			b[i] = M[j+l*M.Dim(COL)];

	int  ldb    = M.Dim(LIN);

	T    beta   = T(0);


	res.Dim(0)  = M.Dim(COL);
	res.Dim(1)  = _dim[LIN];
	res.Reset();
	

	T*   c      = new T[res.Size()];
	int  ldc    = _dim[LIN];

	if (typeid(T) == typeid(double))
		dgemm_ (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
	else if (typeid(T) == typeid(raw))
		cgemm_ (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);

	i = 0;
	for (j = 0; j < res.Dim(COL); j++)
		for (l = 0; l < res.Dim(LIN); l++, i++)
			res[j+l*res.Dim(COL)] = c[i];

	delete [] a;
	delete [] b;
	delete [] c;
	
	return res;
	
}



template<class T>
T
Matrix<T>::norm () const {

	T   res   = (T)0.0;

	int n    = Size();
	int incx = 1;

	if (typeid(T)      == typeid(raw)) {
		res = cblas_scnrm2 (n, _M, incx);
	} else if (typeid(T) == typeid(double))
		res = cblas_scnrm2 (n, _M, incx);

	return pow(res,2);

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
