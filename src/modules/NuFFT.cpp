#include "NuFFT.hpp"

using namespace RRStrategy;

std::string sides[3] = {"Nx", "Ny", "Nz"};


NuFFT::NuFFT () {

}


NuFFT::~NuFFT () {

	this->Finalise();

}


RRSModule::error_code 
NuFFT::Finalise () {

	if (m_initialised)
	  nfft::finalize (&m_fplan, &m_iplan);

	delete [] m_N;
	delete [] m_n;

}



RRSModule::error_code
NuFFT::Init () {

	RRSModule::error_code error = OK; 
	m_initialised               = false;

	m_N   = new int[3];
	m_n   = new int[3];

	for (int i = 0; i < 3; i++) {
		m_N[i] = 0; 
		m_n[i] = 0;
	}

	// Dimensions ---------------------------

	Attribute("dim",       &m_dim);

	for (int i = 0; i < m_dim; i++)
		Attribute (sides[i].c_str(),       &m_N[i]);

	Attribute("M",         &m_M);
	Attribute("shots",     &m_shots);

	// --------------------------------------

	Attribute("maxit",   &m_maxit);
	Attribute("epsilon", &m_epsilon);

	// Oversampling -------------------------

	int      m           = 1;
	double   alpha       = 1.0;

	Attribute ("m",       &m);
	Attribute ("alpha",   &alpha);

	for (int i = 0; i < m_dim; i++)
		m_n[i] = ceil (m_N[i]*alpha);

	// --------------------------------------

	printf ("  intialising nfft::init (%i, {%i, %i, %i}, %i, {%i, %i, %i}, %i, *, *, %.9f)\n", 
			m_dim, 
			m_N[0], m_N[1], m_N[2],
			m_M * m_shots,
			m_n[0], m_n[1], m_n[2],
			m,
			m_epsilon);

	nfft::init (m_dim, m_N, m_M*m_shots, m_n, m, &m_fplan, &m_iplan, m_epsilon);

	m_initialised = true;

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
		m_iplan.y[i][0] = m_raw[i].real();
		m_iplan.y[i][1] = m_raw[i].imag();
	}

	// Copy k-space and weights to allocated memory
	memcpy (&(m_fplan.x[0]), &m_kspace[0], m_kspace.Size()*sizeof(double));
	memcpy (&(m_iplan.w[0]), &m_helper[0], m_helper.Size()*sizeof(double));

	// Precompute PSI
	nfft::weights (&m_fplan, &m_iplan);

	// Inverse FT
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

