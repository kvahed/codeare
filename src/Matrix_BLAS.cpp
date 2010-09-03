template<class T>
Matrix<T> Matrix<T>::dot (Matrix<T>& M) {
	
	Matrix<T> res;

	char transa = 'N';
	char transb = 'N';

	int i = 0, j = 0, l = 0;

	int  m      = 3;
	int  n      = 8;

	int  k      = 8;//_dim[1];

	T    alpha  = T(1);

	T*   a      = new T[Size()];
	i = 0;
	for (j = 0; j < _dim[COL]; j++)
		for (l = 0; l < _dim[LIN]; l++, i++) {
			a[i] = _M[j+l*_dim[COL]];
		}

	int  lda    = 3;//_dim[0];

	T*   b      = new T[M.Size()];
	i = 0;
	for (j = 0; j < M.Dim(COL); j++)
		for (l = 0; l < M.Dim(LIN); l++, i++) {
			b[i] = M[j+l*M.Dim(COL)];
			std::cout << b[i] << std::endl;
		}

	int  ldb    = 8;//M.Dim(0);

	T    beta   = T(0);


	res.Dim(0)  = 3;//_dim[0];
	res.Dim(1)  = 8;//M.Dim(1);
	res.Reset();
	

	T*   c      = new T[res.Size()];
	int  ldc    = res.Dim(1);

	if (typeid(T) == typeid(double))
		dgemm_ (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
	else if (typeid(T) == typeid(raw))
		cgemm_ (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);

	i = 0;
	for (i = 0; i < res.Size(); i++)
		res[i] = c[i];
		
	delete [] a;
	delete [] b;
	delete [] c;
	
	return res;
	
}
