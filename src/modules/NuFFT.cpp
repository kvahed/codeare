#include "NuFFT.hpp"

using namespace RRStrategy;

NuFFT::NuFFT () {

}

NuFFT::~NuFFT () {

	nfft::finalize (&m_fplan, &m_iplan);

	delete [] m_N;
	delete [] m_n;
}

RRSModule::error_code
NuFFT::Init () {

	RRSModule::error_code error = OK; 

	Attribute("dim",     &m_dim);

	m_N = new int[m_dim];
	m_n = new int[m_dim];

	Attribute("Nx",      &m_N[0]);
	Attribute("Ny",      &m_N[1]);
	Attribute("M",       &m_M);

	Attribute("maxit",   &m_maxit);
	Attribute("epsilon", &m_epsilon);

	int      m           = 6;
	double   alpha       = 2.0;

   	m_n[0] = ceil (m_N[0]*alpha);
	m_n[1] = ceil (m_N[1]*alpha);

	// Number of samples
	m_M = m_helper.Size();

	nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan, &m_iplan, m_epsilon);

	return error;

}

RRSModule::error_code
NuFFT::Process () {

	// Some variables
	RRSModule::error_code error = OK;

	printf ("Processing NuFFT ...\n");
	ticks start = getticks();

	// Copy data from incoming matrix to the nufft input array
	for (int i = 0; i < m_raw.Size(); i++) {
		(m_iplan.y[i])[0] = (m_raw[i]).real();
		(m_iplan.y[i])[1] = (m_raw[i]).imag();
	}

	// Kspace adjustment (Don't know yet why necessary)
	m_kspace = m_kspace / (1/(GAMMA/128*m_N[0]));
	
	// Copy k-space and weights to allocated memory
	memcpy (&(m_fplan.x[0]), &m_kspace[0], m_kspace.Size()*sizeof(double));
	memcpy (&(m_iplan.w[0]), &m_helper[0], m_helper.Size()*sizeof(double));

	// Precompute PSI
	nfft::weights (&m_fplan, &m_iplan);

	nfft::ift     (&m_fplan, &m_iplan, m_maxit, m_epsilon);

	// Resize m_raw for output
	for (int i = 0; i < INVALID_DIM; i++)
		m_raw.Dim(i) = 1;
	for (int i = 0; i < m_dim; i++)
		m_raw.Dim(i) = m_N[i];
	m_raw.Reset();

	// Copy back reconstructed image to outgoing matrix
	for (int i = 0; i < m_raw.Size(); i++)
		m_raw[i] = raw (m_iplan.f_hat_iter[i][0], m_iplan.f_hat_iter[i][1]);

	printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), start) / ClockRate());
	
	return error;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new NuFFT;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

