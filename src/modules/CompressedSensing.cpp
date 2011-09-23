#include "CompressedSensing.hpp"
#include "FFT.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_wavelet2d.h>

using namespace RRStrategy;


RRSModule::error_code
CompressedSensing::Init () {

	printf ("Intialising CompressedSensing ...\n");

	for (size_t i = 0; i < 3; i++)
		m_N[i] = 1;

	int wli = 0;
	int m_fft = 0;

	Attribute ("dim",     &m_dim);
	Attribute ("Nx",      &m_N[0]);
	Attribute ("Ny",      &m_N[1]);
	Attribute ("Nz",      &m_N[2]);
    printf ("  Geometry: %iD (%i,%i,%i)\n", m_dim, m_N[0], m_N[1], m_N[2]);
	
	Attribute ("tvw",     &m_cgparam.tvw);
	Attribute ("xfmw",    &m_cgparam.xfmw);
    printf ("  Weights: TV(%.4f) L1(%.4f)\n", m_cgparam.tvw, m_cgparam.xfmw);
	
	Attribute ("fft",     &m_cgparam.fft);
	printf ("  FFT class: ");
	switch (m_fft) 
		{
		case 0:  printf ("%s", "Cartesian"); break;
		case 1:  printf ("%s", "Non-Cartesian"); break;
		default: printf ("%s", "Cartesian"); m_fft = 0; break;
		}
	printf ("\n");

	Attribute ("cgconv", &m_cgparam.cgconv);
	Attribute ("cgiter", &m_cgparam.cgiter);
	printf ("  Maximum %i NLCG iterations or convergence to m_cgconv", m_cgparam.cgiter, m_cgparam.cgconv);	

	

	m_initialised = true;
	printf ("... done.\n\n");

	return RRSModule::OK;

}


RRSModule::error_code
CompressedSensing::Process () {

	Matrix<cplx>*   data  = m_cplx["data"];
	Matrix<double>* pdf   = m_real["pdf"];

	Matrix<cplx>*   im_dc;
	AddCplx ("im_dc", im_dc  = new Matrix<cplx> (data->Dim()));

	ticks tic = getticks();

	printf ("Fourier transforming ... "); fflush(stdout);
	
	for (int i = 0; i < data->Size(); i++)
		data->At(i) /= pdf->At(i);
	
	(*data) = FFT::Backward(*data);

	printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

	tic = getticks();

	cplx ma = data->Maxabs();
	(*data)  /= ma;

	printf ("Wavelet transforming ... "); fflush(stdout);

	gsl_wavelet           *w    = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
	gsl_wavelet_workspace *work = gsl_wavelet_workspace_alloc (data->Size());

	double* re = (double*) malloc (data->Size()*sizeof(double));
	double* im = (double*) malloc (data->Size()*sizeof(double));
	size_t  l  = data->Dim(0); // Needs to be power of 2 and quadratic

	for (size_t j = 0; j < l; j++)
		for (size_t i = 0; i < l; i++) {
			re [i*l+j] = data->At(j*l+i).real();
			im [i*l+j] = data->At(j*l+i).imag();
		}

	if (!(gsl_wavelet2d_nstransform_forward (w, re, l, l, l, work) == GSL_SUCCESS))
		printf ("Wavelet transfor for real part failed\n.");

	if (!(gsl_wavelet2d_nstransform_forward (w, im, l, l, l, work) == GSL_SUCCESS))
		printf ("Wavelet transfor for imaginary part failed\n.");

	for (size_t j = 0; j < l; j++)
		for (size_t i = 0; i < l; i++) {
			im_dc->At(j*l+i) = cplx(re[i*l+j],im [i*l+j]);
		}

	free (re);
	free (im);
	gsl_wavelet_workspace_free (work);
	
	printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
	tic = getticks();

	for (int i = 0; i < m_csiter; i++) {
		
		NLCG (im_dc, &m_cgparam);
		
		// 	res = fnlCg(res,param);
		// 	im_res = XFM'*res;
		// 	figure(100), imshow(abs(im_res),[]), drawnow

	}
	
	return OK;

}





// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CompressedSensing;
}


extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}


