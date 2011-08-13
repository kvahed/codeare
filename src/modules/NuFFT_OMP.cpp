#include "NuFFT_OMP.hpp"
#include "OMP.hpp"

using namespace RRStrategy;

std::string sides[3] = {"Nx", "Ny", "Nz"};

NuFFT_OMP::NuFFT_OMP () {}


NuFFT_OMP::~NuFFT_OMP () {
	
	this->Finalise();
	
}


RRSModule::error_code
NuFFT_OMP::Finalise () {

	if(m_initialised)	
		for (int i = 0; i < NTHREADS; i++)
			nfft::finalize (&m_fplan[i], &m_iplan[i]);

	delete [] m_N;
	delete [] m_n;

	return OK;

}

RRSModule::error_code
NuFFT_OMP::Init () {

	RRSModule::error_code error = OK; 

	m_initialised = false;

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


	Attribute("M",       &m_M);
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

	// Initialise FT plans ------------------
	
	for (int i = 0; i < NTHREADS; i++)
		nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan[i], &m_iplan[i], m_epsilon);

	// --------------------------------------

	m_initialised = true;

	return error;

}

RRSModule::error_code
NuFFT_OMP::Process () {

	// Some variables
	RRSModule::error_code error = OK;

	printf ("Processing NuFFT_OMP ...\n");
	ticks start = getticks();

	int imgsize = 1;
	for (int i = 0; i < m_dim; i++)
		imgsize *= m_N[i];

	Matrix <raw> tmp;

	for (int i = 0; i < m_dim; i++)
		tmp.Dim(i) = m_N[i];

	// Store all arms separately?
	if (m_verbose)
		tmp.Dim(m_dim) = m_shots + 1;

	tmp.Reset();

#pragma omp parallel default (shared) 
	{
		
		omp_set_num_threads(NTHREADS);
		int tid      = omp_get_thread_num();
		
#pragma omp for

		for (int j = 0; j < m_shots; j++) {
			
			int     os         = j * m_M;

			// Copy data from incoming matrix to the nufft input array
			for (int i = 0; i < m_M; i++) {
				(m_iplan[tid].y[i])[0] = (m_raw[i + os]).real();
				(m_iplan[tid].y[i])[1] = (m_raw[i + os]).imag();
			}

			// Copy k-space and weights to allocated memory
			memcpy (&(m_fplan[tid].x[0]), &m_kspace[os * m_dim], m_M * m_dim * sizeof(double));
			memcpy (&(m_iplan[tid].w[0]), &m_helper[os]        , m_M *         sizeof(double));

			// Precompute PSI
			nfft::weights (&m_fplan[tid], &m_iplan[tid]);
			
			nfft::ift     (&m_fplan[tid], &m_iplan[tid], m_maxit, m_epsilon);

			if (m_verbose)
				for (int i = 0; i < imgsize; i++) {
   				tmp[(j+1) * imgsize + i] = raw(m_iplan[tid].f_hat_iter[i][0], m_iplan[tid].f_hat_iter[i][1]);
					tmp[i] += tmp[(j+1) * imgsize + i];
				}
			
			
		}

	}

	printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), start) / Toolbox::Instance()->ClockRate());

	m_raw = tmp;
	return error;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new NuFFT_OMP;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

