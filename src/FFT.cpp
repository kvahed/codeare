#include "FFT.hpp"
#include <fftw3.h>

Matrix<cplx> 
FFT::Forward (Matrix<cplx>& m)  {
	
	assert (m.Is1D() || m.Is2D() || m.Is3D());
	
    Matrix<cplx> res = m;
	fftwf_plan   p;

	Shift(res);
	
	if (m.Is1D())
		p = fftwf_plan_dft_1d (m.Dim(0),                     (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (m.Is2D())
		p = fftwf_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (m.Is3D())
		p = fftwf_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftwf_execute(p);
	fftwf_destroy_plan(p);

	Shift(res);

    return res / res.Size();
	
}

Matrix<cplx>
FFT::Backward (Matrix<cplx>& m) {

	assert (m.Is1D() || m.Is2D() || m.Is3D());
	
    Matrix<cplx> res = m;
	fftwf_plan   p; 
	
	res = Shift(res);

	if (m.Is1D())
		p = fftwf_plan_dft_1d (m.Dim(0),                     (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (m.Is2D())
		p = fftwf_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (m.Is3D())
		p = fftwf_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftwf_execute(p);
	fftwf_destroy_plan(p);

	res = Shift(res);

	return res;
	
}


Matrix<cplx>
FFT::Shift (Matrix<cplx>& m) {
	
	assert (m.Is1D() || m.Is2D() || m.Is3D());
	
	Matrix<cplx> res  = m;
	
	for (size_t s = 0; s < m.Dim(2); s++)
		for (size_t l = 0; l < m.Dim(1); l++)
			for (size_t c = 0; c < m.Dim(0); c++)
				res.At (c,l,s) *= (float) pow ((float)-1.0, (float)(s+l+c));
	
	return res;

}
