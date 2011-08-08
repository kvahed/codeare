#include <fftw3.h>

template <class T>
Matrix<T> Matrix<T>::fft() const {

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
Matrix<T> Matrix<T>::ifft() const {

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
Matrix<T> Matrix<T>::fftshift (const int d) const {

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
Matrix<T> Matrix<T>::ifftshift (const int d) const {

	return fftshift(d);

}
