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
		m_N[i] = 1; 
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

	// Verbosity ----------------------------

	Attribute ("verbose",   &m_verbose);
	printf ("  verbose feedback: %i \n", m_verbose);

	// Oversampling -------------------------

	int      m           = 1;
	double   alpha       = 1.0;

	Attribute ("m",       &m);
	Attribute ("alpha",   &alpha);

	for (int i = 0; i < m_dim; i++)
		m_n[i] = ceil (m_N[i]*alpha);

	// Initialise FT plans ------------------
	
	for (int i = 0; i < NTHREADS; i++)
		nfft::init (m_dim, m_N, m_M*m_shots/NTHREADS, m_n, m, &m_fplan[i], &m_iplan[i], m_epsilon);

	// --------------------------------------

	printf ("  intialising nfft::init (%i, {%i, %i, %i}, %i, {%i, %i, %i}, %i, *, *, %.9f)\n", 
			m_dim, 
			m_N[0], m_N[1], m_N[2],
			m_M * m_shots/NTHREADS,
			m_n[0], m_n[1], m_n[2],
			m,
			m_epsilon);

	m_initialised = true;

	return error;

}

RRSModule::error_code
NuFFT_OMP::Process () {

	// Some variables
	RRSModule::error_code error = OK;

	printf ("Processing NuFFT_OMP ...\n");
	ticks start = getticks();

	Matrix<cplx>*   data    = m_cplx["data"];
	Matrix<double>* kspace  = m_real["kspace"];
	Matrix<double>* weights = m_real["weights"];

	int imgsize = m_N[0] * m_N[1] * m_N[2];

	Matrix <cplx> tmp;

	for (int i = 0; i < m_dim; i++)
		tmp.Dim(i) = m_N[i];

	// Store all arms separately?
	if (m_verbose)
		tmp.Dim(m_dim) = NTHREADS + 1;

	tmp.Reset();

#pragma omp parallel default (shared) 
	{
		
		omp_set_num_threads(NTHREADS);
		int tid      = omp_get_thread_num();
		
#pragma omp for

		for (int j = 0; j < NTHREADS; j++) {
			
			int     os         = j * m_M * m_shots / NTHREADS;

			// Copy data from incoming matrix to the nufft input array
			for (int i = 0; i < m_M*m_shots/NTHREADS; i++) {
				(m_iplan[tid].y[i])[0] = (data->At(i + os)).real();
				(m_iplan[tid].y[i])[1] = (data->At(i + os)).imag();
			}

			// Copy k-space and weights to allocated memory
			memcpy (&(m_fplan[tid].x[0]),  &kspace->At(os * m_dim), m_dim * m_M * m_shots / NTHREADS * sizeof(double));
			memcpy (&(m_iplan[tid].w[0]), &weights->At(         0),         m_M * m_shots / NTHREADS * sizeof(double));

			// Precompute PSI & IFT
			nfft::weights (&m_fplan[tid], &m_iplan[tid]);
			nfft::ift     (&m_fplan[tid], &m_iplan[tid], m_maxit, m_epsilon);

		}


#pragma omp for schedule (dynamic, imgsize/NTHREADS)

	for (int i = 0; i < imgsize; i++)
		for (int j = 0; j < NTHREADS; j++) {
			if (m_verbose)
				tmp[(j+1) * imgsize + i] = cplx(m_iplan[j].f_hat_iter[i][0], m_iplan[j].f_hat_iter[i][1]);
			tmp[i] += cplx(m_iplan[j].f_hat_iter[i][0], m_iplan[j].f_hat_iter[i][1]);
		}
	


	}
			
	printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), start) / Toolbox::Instance()->ClockRate());

	(*data) = tmp;
	return error;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new NuFFT_OMP;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

