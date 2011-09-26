#include "CompressedSensing.hpp"

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

	printf ("Processing CompressedSensing ...\n");
	ticks csstart = getticks();

	Matrix<cplx>*   data  = m_cplx["data"];
	Matrix<double>* pdf   = m_real["pdf"];

	Matrix<cplx>*   im_dc;
	AddCplx ("im_dc", im_dc  = new Matrix<cplx> (data->Dim()));

	ticks tic = getticks();

	printf ("  Fourier transforming ... "); fflush(stdout);
	
	for (int i = 0; i < data->Size(); i++)
		data->At(i) /= pdf->At(i);
	
	(*data) = FFT::Backward(*data);

	printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

	tic = getticks();

	cplx ma = data->Maxabs();
	(*data)  /= ma;

	printf ("  Wavelet transforming ... "); fflush(stdout);

	(*im_dc) = DWT::Forward (*data);

	printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

	tic = getticks();

	printf ("  Running CG ... "); fflush(stdout);

	for (int i = 0; i < m_csiter; i++) {
		
		NLCG ((*im_dc), (*data), m_cgparam);
		
	}
	
	printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

	printf ("... done. WTime: %.4f seconds.\n\n", elapsed(getticks(), csstart) / Toolbox::Instance()->ClockRate());

	return OK;

}





// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CompressedSensing;
}


extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}


