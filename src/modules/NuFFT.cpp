#include "NuFFT.hpp"

using namespace RRStrategy;

NuFFT::NuFFT () {

}

NuFFT::~NuFFT () {

	free (m_ftk);
	free (m_ftw);
	free (m_ftin);
	free (m_ftout);

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

	m_ftin     = (double*) malloc (2          * m_M    * sizeof(double)); 
	m_ftk      = (double*) malloc (    m_dim  * m_M    * sizeof(double)); 
	m_ftw      = (double*) malloc (             m_M    * sizeof(double)); 
	m_ftout    = (double*) malloc (2 * m_N[0] * m_N[1] * sizeof(double)); 

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
		m_ftin[2*i  ] = (m_raw[i]).real();
		m_ftin[2*i+1] = (m_raw[i]).imag();
	}

	// Kspace adjustment (Don't know yet why necessary)
	m_kspace = m_kspace / (1/(GAMMA/128*m_N[0]));
	
	// Copy k-space and weights to allocated memory
	memcpy (m_ftk, &m_kspace[0], m_kspace.Size()*sizeof(double));
	memcpy (m_ftw, &m_helper[0], m_helper.Size()*sizeof(double));

	// Assign k-space and weights to FT
	nfft::kspace  (&m_fplan,           m_ftk);
	nfft::weights (&m_fplan, &m_iplan, m_ftw);
	nfft::ift     (&m_fplan, &m_iplan, m_ftin, m_ftout, m_maxit, m_epsilon);

	// Resize m_raw for output
	for (int i = 0; i < INVALID_DIM; i++)
		m_raw.Dim(i) = 1;
	for (int i = 0; i < m_dim; i++)
		m_raw.Dim(i) = m_N[i];
	m_raw.Reset();

	// Copy back reconstructed image to outgoing matrix
	for (int i = 0; i < m_raw.Size(); i++)
		m_raw[i] = raw(m_ftout[2*i], m_ftout[2*i+1]); 

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

