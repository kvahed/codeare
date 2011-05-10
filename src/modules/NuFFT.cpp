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

	nfft::init (m_dim, m_N, m_helper.Size(), m_n, m, &m_fplan, &m_iplan, m_epsilon);

	return error;

}

RRSModule::error_code
NuFFT::Process () {

	m_M = m_kspace.Size() / m_dim;

	RRSModule::error_code error = OK;

	Matrix<raw> in = m_raw;

	/*for (int i = 0; i < m_dim; i++)
		m_raw.Dim(i) = m_N[i];

	for (int i = m_dim; i < INVALID_DIM; i++)
	m_raw.Dim(i) = 1;*/

	for (int i = 0; i<INVALID_DIM; i++)
		m_raw.Dim(i) = 1;
	
	m_raw.Dim(0) = m_N[0];
	m_raw.Dim(1) = m_N[1];

	m_raw.Reset();

	double*     ftin     = (double*) malloc (2 *    in.Size() * sizeof(double)); 
	double*     ftk      = (double*) malloc ( m_kspace.Size() * sizeof(double)); 
	double*     ftw      = (double*) malloc ( m_helper.Size() * sizeof(double)); 
	double*     ftout    = (double*) malloc (2 * m_raw.Size() * sizeof(double)); 

	for (int i = 0; i < in.Size(); i++) {
		ftin[2*i  ] = (in[i]).real();
		ftin[2*i+1] = (in[i]).imag();
	}

	m_kspace = m_kspace / (1/(GAMMA/128*m_N[0]));
	
	memcpy (ftk, &m_kspace[0], m_kspace.Size()*sizeof(double));
	memcpy (ftw, &m_helper[0], m_helper.Size()*sizeof(double));

	nfft::kspace  (&m_fplan,           ftk);
	nfft::weights (&m_fplan, &m_iplan, ftw);
	nfft::ift     (&m_fplan, &m_iplan, ftin, ftout, m_maxit, m_epsilon);

	for (int i = 0; i < m_raw.Size(); i++)
		m_raw[i] = raw(ftout[2*i], ftout[2*i+1]); 

	delete [] ftout;
	delete [] ftk;
	delete [] ftin;

	return error;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new NuFFT;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

