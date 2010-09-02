template <class T>
int  Matrix<T>::EIG (Matrix<raw>& ev) {
	
#ifdef HAVE_LAPACK

	char    jobvl = 'N';
	char    jobvr = 'N';
	
	int     n     = _dim[COL];

	T*      a     = new T[Size()];

	for (int i = 0; i < Size(); i++)
		a[i] = _M[i];

	int     lda   = n;
	
	T*       w    = new T[n];

	T*       wi;
	if (typeid(T) == typeid(double))
		wi        = new T[n];
	else 
		wi        = new T[1];

	T*       vl   = new T[1];
	T*       vr   = new T[1];

    int     ldvr  =  1;
    int     ldvl  =  1;
	
	T*      work  = new T[1];

	int     info  =  0;
	int     lwork = -1;

	float*  rwork;
	if (typeid(T) == typeid(raw))
		rwork = new float[2*n];
	else 
		rwork = new float[1];

	if (typeid(T) == typeid(raw))
		cgeev_ (&jobvl, &jobvr, &n, a, &lda, w,     vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
	else if (typeid(T) == typeid(double))
		dgeev_ (&jobvl, &jobvr, &n, a, &lda, w, wi, vl, &ldvl, vr, &ldvr, work, &lwork,        &info);
	
	lwork = (int) (raw(work[0]).real());

	std::cout << "xgeev_ (lwork = " << lwork << ")" << std::endl;

	delete [] work;
	work = new T[lwork];

	if (typeid(T) == typeid(raw))
		cgeev_ (&jobvl, &jobvr, &n, a, &lda, w,     vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
	else if (typeid(T) == typeid(double))
		dgeev_ (&jobvl, &jobvr, &n, a, &lda, w, wi, vl, &ldvl, vr, &ldvr, work, &lwork,        &info);

	std::cout << "xgeev_ (info = " << info << ")" << std::endl;

	ev.Dim(COL) = n;
	ev.Reset();
	
	for (int j = 0; j < n; j++)
		ev[j] = w[j];

	std::cout << ev ;

	delete [] a;
	delete [] w;
	delete [] wi;
	delete [] vl;
	delete [] vr;
	delete [] work;
	delete [] rwork;

	return info;

#else

	return -1;

#endif

}


template<class T>
int Matrix<T>::SVD (Matrix<T>& lsv, Matrix<T>& rsv, Matrix<double>& sv) {

#ifdef HAVE_LAPACK

	char    jobz = 'A';
	
	int     m    = _dim[COL];
	int     n    = _dim[LIN];

	int     lda  = _dim[COL];
	int     ldu  = (m >= n) ? m : n;
	int     ldvt = (m >= n) ? n : m;

	T*      a    = new T[Size()];
	
	for (int i = 0; i < Size(); i++)
		a[i] = _M[i];
	
	float*  sf = new float[1];
	double* sd = new double[1];
	
	if (typeid(T) == typeid(raw)) {
		delete [] sf;
		sf = new float [MIN(_dim[0],_dim[1])];
	} else if (typeid(T) == typeid(double)) {
		delete [] sd;
		sd = new double[MIN(_dim[0],_dim[1])];
	}

	T*      u    = new T[ldu*ldu];
	T*      vt   = new T[ldvt*ldvt];
	T*      work = new T[1];
	
	float* rwork = new float[1];
	
	if (typeid(T) == typeid(raw)) {
		delete [] rwork;
		rwork = new float[MIN(m,n)*MAX(5*MIN(m,n)+7,2*MAX(m,n)+2*MIN(m,n)+1)];
	}

	int    lwork = -1;
	int*   iwork = new int[8 * MIN(_dim[0],_dim[1])];
	int     info = 0;

	if (typeid(T) == typeid(raw))
		cgesdd_ (&jobz, &m, &n, a, &lda, sf, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info);
	else if (typeid(T) == typeid(double))
		dgesdd_ (&jobz, &m, &n, a, &lda, sd, u, &ldu, vt, &ldvt, work, &lwork,        iwork, &info);
	
	lwork = (int) raw(work[0]).real();
	
	std::cout << "xgesdd_ (lwork = " << lwork << ")" << std::endl;
	
	delete [] work;
	work = new T[lwork];
	
	if (typeid(T) == typeid(raw))
		cgesdd_ (&jobz, &m, &n, a, &lda, sf, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info);
	else if (typeid(T) == typeid(double))
		dgesdd_ (&jobz, &m, &n, a, &lda, sd, u, &ldu, vt, &ldvt, work, &lwork,        iwork, &info);
	
	lsv.Dim(0) = ldu;
	lsv.Dim(1) = ldu;
	lsv.Reset();
	for (int j = 0; j < ldu*ldu; j++)
		lsv[j] = u[j];
	
	rsv.Dim(0) = ldvt;
	rsv.Dim(1) = ldvt;
	rsv.Reset();
	for (int k = 0; k < ldvt*ldvt; k++)
		rsv[k] = vt[k];
	
	sv.Dim(0) = MIN(_dim[0],_dim[1]);
	sv.Reset();
	for (int l = 0; l < MIN(_dim[0],_dim[1]); l++)
		sv[l] = (typeid(T) == typeid(raw)) ? sf[l] : sd[l];
	
	std::cout << sv;
	
	delete [] a;
	delete [] sf;
	delete [] sd;
	delete [] u;
	delete [] vt;
	delete [] work;
	delete [] rwork;
	delete [] iwork;
	
	return info;

#else //  HAVE_LAPACK
	
	return -1;

#endif // HAVE_LAPACK
	
} 

