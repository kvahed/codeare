template <class T> inline 
Matrix<T>::Matrix () {

    for (size_t i = 0; i < INVALID_DIM; i++) {
        _dim [i] = 1;
        _res [i] = 1.0;
	}

	_M.resize(Size());

	for (size_t i = 0; i < Size(); i++)
		_M[i] = T(0.0);


}



template <class T> inline 
Matrix<T>::Matrix (const size_t n) {

	_dim [COL] = n;
	_dim [LIN] = n;

    for (size_t i = 2; i < INVALID_DIM; i++)
        _dim [i] = 1;

	for (size_t i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;
	
	_M.resize(n*n);
	
	for (size_t i = 0; i < Size(); i++)
		_M[i] = T(0.0);

}



template <class T> inline 
Matrix<T>::Matrix (const size_t m, const size_t n) {

	_dim [0] = m;
	_dim [1] = n;

    for (size_t i = 2; i < INVALID_DIM; i++) 
        _dim [i] = 1;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;

	_M.resize(m*n);

	for (size_t i = 0; i < Size(); i++)
		_M[0] = T(0.0);

}



template <class T> inline 
Matrix<T>::Matrix (const size_t m, const size_t n, const size_t k) {

	_dim [0] = m;
	_dim [1] = n;
	_dim [2] = k;
	
    for (size_t i = 3; i < INVALID_DIM; i++)
        _dim [i] = 1;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;
	
	_M.resize(m*n*k);
	
	for (size_t i = 0; i < Size(); i++)
		_M[0] = T(0.0);
	
}



template <class T> inline 
Matrix<T>::Matrix (const size_t col, const size_t lin, const size_t cha, const size_t set, 
                   const size_t eco, const size_t phs, const size_t rep, const size_t seg, 
                   const size_t par, const size_t slc, const size_t ida, const size_t idb, 
                   const size_t idc, const size_t idd, const size_t ide, const size_t ave) {
	
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
	
    for (size_t i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;
	
	_M.resize(Size());
	
	for (size_t i = 0; i < Size(); i++)
		_M[0] = T(0.0);

}



template <class T> inline 
Matrix<T>::Matrix (const size_t* dim) {
	
	for (size_t i = 0; i < INVALID_DIM; i++) {
		_dim[i] = dim[i];
        _res[i] = 1.0;
	}
	
	_M.resize(Size());
	
	for (size_t i = 0; i < Size(); i++)
		_M[0] = T(0.0);
	
}



template <class T> inline 
Matrix<T>::Matrix (const Matrix<T> &M) {
	
	for (size_t i = 0; i < INVALID_DIM; i++) {
		_dim[i] = M.Dim(i);
		_res[i] = M.Res(i);
	}
	   
	_M.resize(Size());
	
	memcpy (&_M[0], M.Data(), Size() * sizeof(T));
	
}



template <class T> inline 
Matrix<T>::~Matrix() {
    
#ifdef PARC_MODULE_NAME
    ICE_SET_FN ("Matrix<T>::~Matrix()")
    ICE_WARN   ("Freeing " << (float)Size() * sizeof(T) / 1024 << " kB of RAM.");
#endif

	_M.clear();
    
}



template <class T> inline 
Matrix<T> Matrix<T>::Id (const size_t n) {

 	Matrix<T> M (n);

 	for (size_t i = 0; i < n; i++)
 		M[i*n+i] = T(1.0);

 	return M;

}



template <class T> inline 
Matrix<T> Matrix<T>::Ones (const size_t m, const size_t n, const size_t l) {

 	Matrix<T> M (m,n,l);

 	for (size_t i = 0; i < M.Size(); i++)
 		M[i] = T(1.0);

 	return M;

}



template <class T> inline 
Matrix<T> Matrix<T>::Ones (const size_t m, const size_t n) {

 	Matrix<T> M (m,n);

 	for (size_t i = 0; i < M.Size(); i++)
 		M[i] = T(1.0);

 	return M;

}



template <class T> inline 
Matrix<T> Matrix<T>::Ones (const size_t n) {

 	return Ones(n,n);

}



template <class T> inline 
Matrix<T> Matrix<T>::Zeros (const size_t n, const size_t m, const size_t l) {

 	Matrix<T> M (m,n,l);

 	for (size_t i = 0; i < M.Size(); i++)
 		M[i] = T(0.0);

 	return M;

}



template <class T> inline 
Matrix<T> Matrix<T>::Zeros (const size_t n, const size_t m) {

 	Matrix<T> M (m,n);

 	for (size_t i = 0; i < M.Size(); i++)
 		M[i] = T(0.0);

 	return M;

}



template <class T> inline 
Matrix<T> Matrix<T>::Zeros (const size_t n) {

 	return Zeros(n,n);

}



template <class T> inline 
Matrix<T> Matrix<T>::Circle (const float* p, const size_t n) {

	Matrix<T> res = Matrix<T>::Zeros(n);

	float m[2];
	float rad;

	rad = p[0] * float(n) / 2.0;

	m[0] = (1.0 - p[1]) * float(n) / 2.0;
	m[1] = (1.0 - p[2]) * float(n) / 2.0;

	for (size_t r = 0; r < res.Dim(1); r++)
		for (size_t c = 0; c < res.Dim(0); c++)
			res(c,r) = ( pow(((float)c-m[0])/rad, 2.0 ) + pow(((float)r-m[0])/rad, 2.0) <= 1.0) ? T(1.0) : T(0.0);

	return res;

}



template <class T> inline 
Matrix<T> Matrix<T>::Sphere (const float* p, const size_t n) {

	Matrix<T> res = Matrix<T>::Zeros(n,n,n);

	float m[3];
	float rad;

	rad = p[0] * float(n) / 2.0;

	m[0] = (1.0 - p[1]) * float(n) / 2.0;
	m[1] = (1.0 - p[2]) * float(n) / 2.0;
	m[2] = (1.0 - p[3]) * float(n) / 2.0;

	for (size_t s = 0; s < res.Dim(2); s++)
		for (size_t r = 0; r < res.Dim(1); r++)
			for (size_t c = 0; c < res.Dim(0); c++)
				res(c,r) = ( pow (((float)c-m[0])/rad, 2.0) + pow (((float)r-m[1])/rad, 2.0) + pow (((float)s-m[2])/rad, 2.0) <= 1.0) ? T(1.0) : T(0.0);

	return res;

}



template <class T> inline 
Matrix<T> Matrix<T>::Ellipse (const float* p, const size_t n, const T v) {

	Matrix<T> res = Matrix<T>::Zeros(n);

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
		
		size_t tid      = omp_get_thread_num();
		size_t chunk    = n / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk) 
		
	for (size_t r = 0; r < n; r++)
		for (size_t c = 0; c < n; c++)
			res(c,r) = (pow( (((float)c-m[1])*cosp+((float)r-m[0])*sinp)/a[1], 2.0 ) + 
						pow( (((float)r-m[0])*cosp-((float)c-m[1])*sinp)/a[0], 2.0) <= 1.0) ? v : T(0.0);

	}

	return res;

}



template <class T> inline 
Matrix<T> Matrix<T>::Ellipsoid (const float* p, const size_t n, const T v) {

	Matrix<T> res = Matrix<T>::Zeros(n,n,n);

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
		
		size_t tid      = omp_get_thread_num();
		size_t chunk    = n / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk) 
		
		for (size_t s = 0; s < n; s++)
			for (size_t r = 0; r < n; r++)
				for (size_t c = 0; c < n; c++)
					res(c,r,s) = ( pow( (((float)c-m[1])*cosp+((float)r-m[0])*sinp)/a[1], 2.0) + 
								   pow( (((float)r-m[0])*cosp-((float)c-m[1])*sinp)/a[0], 2.0) +
								   pow( ((float)s-m[2])/a[2], 2.0) <= 1.0) ? v : T(0.0);
		
	}

	return res;

}



template <class T> inline 
Matrix<T> Matrix<T>::Phantom2D (const size_t n) {

	const size_t ne = 10; // Number of ellipses
	const size_t np = 5;  // Number of geometrical parameters

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

	// Size_Tensities
	T v[ne] = {T(1.0), T(-0.8), T(-0.2), T(-0.2), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1)};

	// Empty matrix
	Matrix<T> res = Matrix<T>::Zeros(n);
	Matrix<T>        e;

	for (size_t i = 0; i < ne; i++) {
		e    = Matrix<T>::Ellipse (p[i], n, v[i]);
		res += e;
	}

	return res;

}



template <class T> inline 
Matrix<T> Matrix<T>::Phantom3D (const size_t n) {

	const size_t ne = 10; // Number of ellipses
	const size_t np =  9; // Number of geometrical parameters

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

	Matrix<T> res = Matrix<T>::Zeros(n,n,n);
	Matrix<T> e;
	
	for (size_t i = 0; i < ne; i++) {
		e    = Matrix<T>::Ellipsoid (p[i], n, v[i]);
		res += e;
	}

	return res;

}


template<> inline Matrix<double> 
Matrix<cplx>::Real () const {

	Matrix<double> res (_dim);

#pragma omp parallel default (shared) 
	{
		
		size_t tid      = omp_get_thread_num();
		size_t chunk    = Size() / omp_get_num_threads();
		
#pragma omp for 
		
		for (size_t i = 0; i < Size(); i++)
			res[i] = _M[i].real();

	}		
		
	return res;

}
    

template<> inline Matrix<double> 
Matrix<cplx>::Imag () const {

	Matrix<double> res (_dim);

#pragma omp parallel default (shared) 
	{
		
		size_t tid      = omp_get_thread_num();
		size_t chunk    = Size() / omp_get_num_threads();
		
#pragma omp for 
		
		for (size_t i = 0; i < Size(); i++)
			res[i] = _M[i].imag();

	}		
		
	return res;

}
    

template<class T> Matrix<size_t>
Matrix<T>::MeshGrid (const Matrix<size_t>& d) {

	size_t side [3];

	side[0] = d(0,1) - d(0,0) + 1;
	side[1] = d(1,1) - d(1,0) + 1;
	side[2] = d(2,1) - d(2,0) + 1;

    Matrix<size_t> mg (side[1], side[0], side[2], 3);

	std::cout << mg.DimsToCString() << std::endl;

	for (size_t s = 0; s < side[2]; s++)
		for (size_t l = 0; l < side[0]; l++)
			for (size_t c = 0; c < side[1]; c++) {
				mg(c,l,s,0) = l + d(0,0);
				mg(c,l,s,1) = c + d(0,1);
				mg(c,l,s,2) = s + d(0,2);
			}
	
	return mg;

}
