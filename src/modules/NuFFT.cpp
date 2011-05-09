#include "NuFFT.hpp"

using namespace RRStrategy;

NuFFT::NuFFT () {

	Attribute("dim",     &m_dim);

	m_N = new int[m_dim];
	m_n = new int[m_dim];

	Attribute("Nx",      &m_N[0]);
	Attribute("Ny",      &m_N[1]);

	Attribute("maxit",   &m_maxit);
	Attribute("epsilon", &m_epsilon);

	int      m           = 6;
	double   alpha       = 2.0;

	m_n[0] = ceil (m_N[0]*alpha);
	m_n[1] = ceil (m_N[1]*alpha);

	nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan, &m_iplan, m_epsilon);

}

NuFFT::~NuFFT () {
	delete [] m_N;
	delete [] m_n;
}

RRSModule::error_code
NuFFT::Process () {

	m_M = m_helper.Size() / m_dim;

	// Check: m_helper.Size() = m_dim * m_raw.Size()
	m_raw.dump (std::string("m_raw.h5").c_str());
	m_helper.dump (std::string("m_helper.h5").c_str());
	DumpConfig (std::string("m_config.xml").c_str());

	RRSModule::error_code error = OK;

	Matrix<raw> in = m_raw;

	for (int i = 0; i < m_dim; i++)
		m_raw.Dim(i) = m_N[i];

	for (int i = m_dim; i < INVALID_DIM; i++)
		m_raw.Dim(i) = 1;
	
	m_raw.Reset();

	double*     ftin     = (double*) malloc (2 *    in.Size() * sizeof(double)); 
	double*     ftk      = (double*) malloc ( m_helper.Size() * sizeof(double)); 
	double*     ftout    = (double*) malloc (2 * m_raw.Size() * sizeof(double)); 

	/*for (int i = 0; i < in.Size(); i++) {
		ftin[2*i  ] = (in[i]).real();
		ftin[2*i+1] = (in[i]).imag();
	}

	memcpy (ftk, &m_helper[0], m_helper.Size()*sizeof(double));
	nfft::kspace (&m_fplan, ftk);

	nfft::ift    (&m_fplan, &m_iplan, ftin, ftout, m_maxit, m_epsilon);


	for (int i = 0; i < m_raw.Size(); i++)
		m_raw[i] = raw(ftout[2*i], ftout[2*i+1]); 
	*/

	m_raw.dump(std::string("m_raw.h5").c_str());

	delete [] ftout;
	delete [] ftk;
	delete [] ftin;

	return error;

}

