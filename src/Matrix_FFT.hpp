#include <fftw3.h>

template <class T>
Matrix<T> Matrix<T>::FFT() const {
	
	assert (Is1D() || Is2D() || Is3D());
	assert (typeid(T) == typeid(raw));
	
    Matrix<T>  res (_dim);
	fftwf_plan p; 
	
#ifdef DEBUG
	if (Is1D()) 
		printf ("Forward 1D FFT: %i\n", _dim[0]);
	else if (Is2D())
		printf ("Forward 2D FFT: %i\n", _dim[0]*_dim[1]);
	else if (Is3D())
		printf ("Forward 3D FFT: %i\n", _dim[0]*_dim[1]*_dim[2]);
#endif
	
	// Data is column-major (inverse order of indices)
	if (Is1D())
		p = fftwf_plan_dft_1d(_dim[0],                   reinterpret_cast<fftwf_complex*>(&_M[0]), reinterpret_cast<fftwf_complex*>(&res[0]), FFTW_FORWARD, FFTW_ESTIMATE);
	else if (Is2D())
		p = fftwf_plan_dft_2d(_dim[1], _dim[0],          reinterpret_cast<fftwf_complex*>(&_M[0]), reinterpret_cast<fftwf_complex*>(&res[0]), FFTW_FORWARD, FFTW_ESTIMATE);
	else if (Is3D())
		p = fftwf_plan_dft_3d(_dim[2], _dim[1], _dim[0], reinterpret_cast<fftwf_complex*>(&_M[0]), reinterpret_cast<fftwf_complex*>(&res[0]), FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftwf_execute(p);
	fftwf_destroy_plan(p);
	
    return res / Size();
	
}

template <class T>
Matrix<T> Matrix<T>::IFFT() const {
	
	assert (Is1D() || Is2D() || Is3D());
	assert (typeid(T) == typeid(raw));
	
    Matrix<T>  res (_dim);
	fftwf_plan p; 
	
#ifdef DEBUG
	if (Is1D()) 
		printf ("Backward 1D FFT: %i\n", _dim[0]);
	else if (Is2D())
		printf ("Backward 2D FFT: %i\n", _dim[0]*_dim[1]);
	else if (Is3D())
		printf ("Backward 3D FFT: %i\n", _dim[0]*_dim[1]*_dim[2]);
#endif
	
	// Data is column-major (inverse order of indices)
	if (Is1D())
		p = fftwf_plan_dft_1d(_dim[0],                   reinterpret_cast<fftwf_complex*>(&_M[0]), reinterpret_cast<fftwf_complex*>(&res[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (Is2D())
		p = fftwf_plan_dft_2d(_dim[1], _dim[0],          reinterpret_cast<fftwf_complex*>(&_M[0]), reinterpret_cast<fftwf_complex*>(&res[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (Is3D())
		p = fftwf_plan_dft_3d(_dim[2], _dim[1], _dim[0], reinterpret_cast<fftwf_complex*>(&_M[0]), reinterpret_cast<fftwf_complex*>(&res[0]), FFTW_BACKWARD, FFTW_ESTIMATE);

	fftwf_execute(p);
	fftwf_destroy_plan(p);
	
    return res;
	
}


template <class T>
Matrix<T> Matrix<T>::FFTShift (const int d) const {
	
	assert (Is1D() || Is2D() || Is3D());
	assert (typeid(T) == typeid(raw));
	
	Matrix<T> res (_dim);
	
	for (int s = 0; s < _dim[2]; s++)
		for (int l = 0; l < _dim[1]; l++)
			for (int c = 0; c < _dim[0]; c++)
				res.At (c,l,s) = this->At(c,l,s) * (float) pow ((float)-1.0, (float)(s+l+c));
	
	return res;
	
}


template <class T>
Matrix<T> Matrix<T>::IFFTShift (const int d) const {
	
	return FFTShift(d);
	
}


template <class T>
Matrix<T> Matrix<T>::HannWindow (const int dim) const {
	
	assert (Is1D() || Is2D() || Is3D());
	
	Matrix<T> res;
	float     h, d;
	float     m[3];
	
	if (Is1D()) {
		
		m[0] = 0.5 * (float)_dim[0];
		m[1] = 0.0;
		m[2] = 0.0;
		
	} else if (Is2D()) {
		
		m[0] = 0.5 * (float)_dim[0];
		m[1] = 0.5 * (float)_dim[1];
		m[2] = 0.0;
		
	} else {

		m[0] = 0.5 * (float)_dim[0];
		m[1] = 0.5 * (float)_dim[1];
		m[2] = 0.5 * (float)_dim[2];

	}
	
	res = this->Squeeze();
	
	for (int s = 0; s < _dim[2]; s++)
		for (int r = 0; r < _dim[1]; r++)
			for (int c = 0; c < _dim[0]; c++) {
				d = pow( pow(((float)c-m[0])/m[0],2.0) + pow(((float)r-m[1])/m[1],2.0) + pow(((float)s-m[2])/m[2],2.0) , 0.5); 
				h = (d < 1) ? (0.5 + 0.5 * cos (PI * d)) : 0.0;
				res(c,r,s) = this->At(c,r,s) * h;
			}
	
	return res;
	
}


template <class T>
Matrix<T> Matrix<T>::SOS (const int d) const {
	
	assert (_dim[d] > 1);
	
	unsigned short nd = this->HDim();
	int dim [INVALID_DIM];
	
	for (int i = 0; i < INVALID_DIM; i++)
		dim[i] = (i != nd) ? _dim[i] : 1;
	
	Matrix<T> res (dim);
	
#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
		for (int i = 0; i < res.Size(); i++) {
			for (int j = 0; j < _dim[nd]; j++)
				res [i] = pow (_M[i + j * res.Size()], 2.0);
			pow (res[i],0.5);
		}

	}

	return res;

}


template <class T> void
Matrix<T>::Squeeze () {
	
	int found = 0;
	
	for (int i = 0; i < INVALID_DIM; i++)
		if (_dim[i] > 1)
			_dim[found++] = _dim[i];
	
	for (int i = found; i < INVALID_DIM; i++)
		_dim[i] = 1;
	
}


template <class T> Matrix<T>
Matrix<T>::Squeeze () const {
	
	Matrix<T> res = (*this);
	res.Squeeze();
	return res;
	
}


template <class T> void
Matrix<T>::Mean (const int d) {

	float quot = (float) d;

	this->Sum (d);
	
	for (int i = 0; i < Size(); i++)
		_M[i] = _M[i] / quot;
	
}


template <class T> Matrix<T>
Matrix<T>::Mean (const int d) const {
	
	Matrix<T> res = (*this);
	res.Mean(d);
	return res;
	
}


template <class T> void
Matrix<T>::Sum (const int d) {

	assert (d>=0 && d < INVALID_DIM);
   
	// No meaningful sum over particular dimension
	if (_dim[d] == 1)
		return;

	// No RAM allocation 
	if (!nb_alloc)
		return;

	// Save old data and resize matrix 
	T* tmp = (T*) malloc (Size() * sizeof (T));
	memcpy (tmp, _M, Size() * sizeof (T));
	free (_M);
	_M       = (T*) malloc ((Size() / _dim[d]) * sizeof(T));

	// Inner size 
	int insize = 1;
	for (int i = 0; i < d; i++)
		insize *= _dim[i];
	
	// Outer size
	int outsize = 1;
	for (int i = d+1; i < INVALID_DIM; i++)
		outsize *= _dim[i];
	
	// Sum
#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		
		for (int i = 0; i < outsize; i++) {
			
#pragma omp for

			for (int j = 0; j < insize; j++) {
				_M[j+i*insize] = 0;
				for (int k = 0; k < _dim[d]; k++)
					_M[j+i*insize] += tmp[j+i*insize*_dim[d]+k*insize];
			}
			
		}

	}
	// Adjust dminesions and clear tmp
	_dim[d] = 1;
	free (tmp);
	
}


template <class T> Matrix<T>
Matrix<T>::Sum (const int d) const {
	
	Matrix<T> res = (*this);
	res.Sum(d);
	return res;
	
}


template <class T> unsigned short
Matrix<T>::HDim () const {
	
	unsigned short nd = 0;
	
	for (int i = 0; i < INVALID_DIM; i++)
		nd  = (_dim[i] > 1) ? i : nd;
	
	return nd;

}


template <class T> void
Matrix<T>::PrintDims () const {

	printf ("Dimensions: ");
	
	for (int i = 0; i < INVALID_DIM; i++)

		printf ("%i ", _dim[i]);
	
	printf ("\n");
	
}

