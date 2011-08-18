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

    for (int i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;

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

    for (int i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;

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

    for (int i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;

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

    for (int i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;

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

    for (int i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;

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
	
    for (int i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;

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
Matrix<T> Matrix<T>::Ones (const int m, const int n, const int l) {

 	static Matrix<T> M (m,n,l);

 	for (int i = 0; i < M.Size(); i++)
 		M[i] = T(1.0);

 	return M;

}



template <class T>
Matrix<T> Matrix<T>::Ones (const int m, const int n) {

 	static Matrix<T> M (m,n);

 	for (int i = 0; i < M.Size(); i++)
 		M[i] = T(1.0);

 	return M;

}



template <class T>
Matrix<T> Matrix<T>::Ones (const int n) {

 	return Ones(n,n);

}



template <class T>
Matrix<T> Matrix<T>::Zeros (const int n, const int m, const int l) {

 	static Matrix<T> M (m,n,l);

 	for (int i = 0; i < M.Size(); i++)
 		M[i] = T(0.0);

 	return M;

}



template <class T>
Matrix<T> Matrix<T>::Zeros (const int n, const int m) {

 	static Matrix<T> M (m,n);

 	for (int i = 0; i < M.Size(); i++)
 		M[i] = T(0.0);

 	return M;

}



template <class T>
Matrix<T> Matrix<T>::Zeros (const int n) {

 	return Zeros(n,n);

}



template <class T>
Matrix<T> Matrix<T>::Circle (const float* p, const int n) {

	static Matrix<T> res = Matrix<T>::Zeros(n);

	float m[2];
	float rad;

	rad = p[0] * float(n) / 2.0;

	m[0] = (1.0 - p[1]) * float(n) / 2.0;
	m[1] = (1.0 - p[2]) * float(n) / 2.0;

	for (int r = 0; r < res.Dim(1); r++)
		for (int c = 0; c < res.Dim(0); c++)
			res(c,r) = ( pow(((float)c-m[0])/rad, 2.0 ) + pow(((float)r-m[0])/rad, 2.0) <= 1.0) ? T(1.0) : T(0.0);

	return res;

}



template <class T>
Matrix<T> Matrix<T>::Sphere (const float* p, const int n) {

	static Matrix<T> res = Matrix<T>::Zeros(n,n,n);

	float m[3];
	float rad;

	rad = p[0] * float(n) / 2.0;

	m[0] = (1.0 - p[1]) * float(n) / 2.0;
	m[1] = (1.0 - p[2]) * float(n) / 2.0;
	m[2] = (1.0 - p[3]) * float(n) / 2.0;

	for (int s = 0; s < res.Dim(2); s++)
		for (int r = 0; r < res.Dim(1); r++)
			for (int c = 0; c < res.Dim(0); c++)
				res(c,r) = ( pow (((float)c-m[0])/rad, 2.0) + pow (((float)r-m[1])/rad, 2.0) + pow (((float)s-m[2])/rad, 2.0) <= 1.0) ? T(1.0) : T(0.0);

	return res;

}



template <class T>
Matrix<T> Matrix<T>::Ellipse (const float* p, const int n, const T v) {

	static Matrix<T> res = Matrix<T>::Zeros(n);

	float m[2];
	float a[2];

	a[0] = p[0] * float(n) / 2.0;
	a[1] = p[1] * float(n) / 2.0;

	m[0] = (1.0 - p[2]) * float(n) / 2.0;
	m[1] = (1.0 - p[3]) * float(n) / 2.0;

	float cosp = cos(p[4]);
	float sinp = sin(p[4]);
	
#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = n / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk) 
		
	for (int r = 0; r < n; r++)
		for (int c = 0; c < n; c++)
			res(c,r) = (pow( (((float)c-m[1])*cosp+((float)r-m[0])*sinp)/a[1], 2.0 ) + 
						pow( (((float)r-m[0])*cosp-((float)c-m[1])*sinp)/a[0], 2.0) <= 1.0) ? v : T(0.0);

	}

	return res;

}



template <class T>
Matrix<T> Matrix<T>::Ellipsoid (const float* p, const int n, const T v) {

	static Matrix<T> res = Matrix<T>::Zeros(n,n,n);

	float m[3];
	float a[3];
	float d;

	a[0] = p[0] * float(n) / 2.0;
	a[1] = p[1] * float(n) / 2.0;
	a[2] = p[2] * float(n) / 2.0;

	m[0] = (1.0 - p[3]) * float(n) / 2.0;
	m[1] = (1.0 - p[4]) * float(n) / 2.0;
	m[2] = (1.0 - p[5]) * float(n) / 2.0;

	float cosp = cos(p[6]);
	float sinp = sin(p[6]);
	
#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = n / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk) 
		
		for (int s = 0; s < n; s++)
			for (int r = 0; r < n; r++)
				for (int c = 0; c < n; c++)
					res(c,r,s) = ( pow( (((float)c-m[1])*cosp+((float)r-m[0])*sinp)/a[1], 2.0) + 
								   pow( (((float)r-m[0])*cosp-((float)c-m[1])*sinp)/a[0], 2.0) +
								   pow( ((float)s-m[2])/a[2], 2.0) <= 1.0) ? v : T(0.0);
		
	}

	return res;

}



template <class T>
Matrix<T> Matrix<T>::Phantom2D (const int n) {

	const int ne = 10; // Number of ellipses
	const int np = 5;  // Number of geometrical parameters

	float p[ne][np] = {
		{ 0.6900, 0.9200,  0.00,  0.0000,  0.0 },
		{ 0.6624, 0.8740,  0.00, -0.0184,  0.0 },
        { 0.1100, 0.3100, -0.22,  0.0000, -0.3 },
		{ 0.1600, 0.4100,  0.22,  0.0000,  0.3 },
		{ 0.2100, 0.2500,  0.00,  0.3500,  0.0 },
		{ 0.0460, 0.0460,  0.00,  0.1000,  0.0 },
		{ 0.0460, 0.0460,  0.00, -0.1000,  0.0 },
		{ 0.0460, 0.0230,  0.08, -0.6050,  0.0 },
		{ 0.0230, 0.0230,  0.00, -0.6060,  0.0 },
		{ 0.0230, 0.0460, -0.06, -0.6050,  0.0 }
	};

	// Intensities
	T v[ne] = {T(1.0), T(-0.8), T(-0.2), T(-0.2), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1)};

	// Empty matrix
	static Matrix<T> res = Matrix<T>::Zeros(n);
	Matrix<T>        e;

	for (int i = 0; i < ne; i++) {
		e    = Matrix<T>::Ellipse (p[i], n, v[i]);
		res += e;
	}

	return res;

}



template <class T>
Matrix<T> Matrix<T>::Phantom3D (const int n) {

	const int ne = 10; // Number of ellipses
	const int np =  9; // Number of geometrical parameters

	float p[ne][np] = {
		{ 0.690, 0.920, 0.900,  0.00,  0.000,  0.000,  0.0, 0.0, 0.0 },
        { 0.662, 0.874, 0.880,  0.00,  0.000,  0.000,  0.0, 0.0, 0.0 },
        { 0.110, 0.310, 0.220, -0.22,  0.000, -0.250, -0.3, 0.0, 0.0 },
        { 0.160, 0.410, 0.210,  0.22,  0.000, -0.250,  0.3, 0.0, 0.0 },
        { 0.210, 0.250, 0.500,  0.00,  0.350, -0.250,  0.0, 0.0, 0.0 },
        { 0.046, 0.046, 0.046,  0.00,  0.100, -0.250,  0.0, 0.0, 0.0 },
        { 0.046, 0.023, 0.020,  0.08, -0.650, -0.250,  0.0, 0.0, 0.0 },
        { 0.046, 0.023, 0.020,  0.06, -0.650, -0.250,  0.0, 0.0, 0.0 },
        { 0.056, 0.040, 0.100, -0.06, -0.105,  0.625,  0.0, 0.0, 0.0 },
        { 0.056, 0.056, 0.100,  0.00,  0.100,  0.625,  0.0, 0.0, 0.0 }
	};

	T v[ne] = {2.0, -0.8, -0.2, -0.2, 0.2, 0.2, 0.1, 0.1, 0.2, -0.2};

	static Matrix<T> res = Matrix<T>::Zeros(n,n,n);
	Matrix<T> e;
	
	for (int i = 0; i < ne; i++) {
		e    = Matrix<T>::Ellipsoid (p[i], n, v[i]);
		res += e;
	}

	return res;

}


/*template<>
Matrix<double> Matrix<raw>::Real () const {

	Matrix<double> res (_dim);

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for 
		
		for (int i = 0; i < Size(); i++)
			res[i] = _M[i].real();

	}		
		
	return res;

	}*/
    
/*template<>
Matrix<double> Matrix<raw>::Imag () const {

	Matrix<double> res (_dim);

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for 
		
		for (int i = 0; i < Size(); i++)
			res[i] = _M[i].imag();

	}		
		
	return res;

	}*/
    
