template <class T> 
Matrix<T>::Matrix () {

	nb_alloc = 0;

    for (int i = 0; i < INVALID_DIM; i++)
        _dim [i] = 1;

    _M = (T*) malloc (Size()*sizeof(T));

	for (int i = 0; i < Size(); i++)
		_M[i] = T(0.0);

    nb_alloc++;


}


template <class T> 
Matrix<T>::Matrix (const int n) {

	nb_alloc = 0;

	_dim [COL] = n;
	_dim [LIN] = n;

    for (int i = CHA; i < INVALID_DIM; i++)
        _dim [i] = 1;

    _M = (T*) malloc (n*n*sizeof(T));

	for (int i = 0; i < Size(); i++)
		_M[i] = T(0.0);

    nb_alloc++;

}


template <class T> 
Matrix<T>::Matrix (const int m, const int n) {

	nb_alloc = 0;

	_dim [0] = m;
	_dim [1] = n;

    for (int i = 2; i < INVALID_DIM; i++)
        _dim [i] = 1;

    _M = (T*) malloc (Size()*sizeof(T));

	for (int i = 0; i < Size(); i++)
		_M[0] = T(0.0);

    nb_alloc++;

}


template <class T> 
Matrix<T>::Matrix (const int m, const int n, const int k) {

	nb_alloc = 0;

	_dim [0] = m;
	_dim [1] = n;
	_dim [2] = k;

    for (int i = 3; i < INVALID_DIM; i++)
        _dim [i] = 1;

    _M = (T*) malloc (Size()*sizeof(T));

	for (int i = 0; i < Size(); i++)
		_M[0] = T(0.0);

    nb_alloc++;

}


template <class T> 
Matrix<T>::Matrix (const int m, const int n, const int k, const int l) {

	nb_alloc = 0;

	_dim [0] = m;
	_dim [1] = n;
	_dim [2] = k;
	_dim [3] = l;
	
    for (int i = 4; i < INVALID_DIM; i++)
        _dim [i] = 1;

    _M = (T*) malloc (Size()*sizeof(T));

	for (int i = 0; i < Size(); i++)
		_M[0] = T(0.0);

    nb_alloc++;

}


template <class T>
Matrix<T>::Matrix (const int col, const int lin, const int cha, const int set, 
                   const int eco, const int phs, const int rep, const int seg, 
                   const int par, const int slc, const int ida, const int idb, 
                   const int idc, const int idd, const int ide, const int ave) {

	nb_alloc = 0;

    _dim[COL] = col;
    _dim[LIN] = lin;
    _dim[CHA] = cha;
    _dim[SET] = set;
    _dim[ECO] = eco;
    _dim[PHS] = phs;
    _dim[REP] = rep;
    _dim[SEG] = seg;
    _dim[PAR] = par;
    _dim[SLC] = slc;
    _dim[IDA] = ida;
    _dim[IDB] = idb;
    _dim[IDC] = idc;
    _dim[IDD] = idd;
    _dim[IDE] = ide;
    _dim[AVE] = ave;

    _M = (T*) malloc (Size() * sizeof (T));

	for (int i = 0; i < Size(); i++)
		_M[0] = T(0.0);

    nb_alloc++;


}


template <class T>
Matrix<T>::Matrix (const int* dim) {

	nb_alloc = 0;

	for (int i = 0; i < INVALID_DIM; i++)
		_dim[i] = dim[i];

    _M = (T*) malloc (Size() * sizeof (T));

	for (int i = 0; i < Size(); i++)
		_M[0] = T(0.0);

    nb_alloc++;

}


template <class T>
Matrix<T>::Matrix (const Matrix<T> &M) {
	
	nb_alloc = 0;
	
	for (int i = 0; i < INVALID_DIM; i++) 
		_dim[i] = M.Dim(i);
	
	_M = (T*) malloc (Size() * sizeof (T));
	
	memcpy (_M, M.Data(), Size() * sizeof(T));
	
	nb_alloc++;
	
}


template <class T> 
Matrix<T>::~Matrix() {
    
#ifdef PARC_MODULE_NAME
    ICE_SET_FN ("Matrix<T>::~Matrix()")
    ICE_WARN   ("Freeing " << (float)Size() * sizeof(T) / 1024 << " kB of RAM.");
#endif

    if (nb_alloc) {
    	free (_M);
        nb_alloc--;
    }
    
}


template <class T>
Matrix<T> Matrix<T>::Id (const int n) {

 	static Matrix<T> M (n);

 	for (int i = 0; i < n; i++)
 		M[i*n+i] = T(1.0);

 	return M;

}


template <class T>
Matrix<T> Matrix<T>::Ones (const int m, const int n) {

 	static Matrix<T> M (m,n);

 	for (int i = 0; i < Size(); i++)
 		M[i] = T(1.0);

 	return M;

}


template <class T>
Matrix<T> Matrix<T>::Ones (const int n) {

 	return Ones(n,n);

}


template <class T>
Matrix<T> Matrix<T>::Zeros (const int n, const int m) {

 	static Matrix<T> M (m,n);

 	for (int i = 0; i < Size(); i++)
 		M[i] = T(0.0);

 	return M;

}


template <class T>
Matrix<T> Matrix<T>::Zeros (const int n) {

 	return Zeros(n,n);

}


