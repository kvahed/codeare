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
Matrix<T> Matrix<T>::HannWindow (const int d) const {
	
	assert (Is1D() || Is2D() || Is3D());
	
	Matrix<T> res (_dim);
	
	//for ()
	
	//x = (0:m-1)'/(n-1);
	//w = 0.5 - 0.5*cos(2*pi*x);
	
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
