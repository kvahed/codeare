template<class T>
Matrix<T> Matrix<T>::dot (Matrix<T>& M) {
	
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

	int  lda    = _dim[LIN];//_dim[0];

	T*   b      = new T[M.Size()];
	i = 0;
	for (j = 0; j < M.Dim(COL); j++)
		for (l = 0; l < M.Dim(LIN); l++, i++)
			b[i] = M[j+l*M.Dim(COL)];

	int  ldb    = M.Dim(LIN);//M.Dim(0);

	T    beta   = T(0);


	res.Dim(0)  = M.Dim(COL);//_dim[0];
	res.Dim(1)  = _dim[LIN];//M.Dim(1);
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
