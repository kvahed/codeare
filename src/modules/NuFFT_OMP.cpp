#include "NuFFT_OMP.hpp"
#include "OMP.hpp"

using namespace RRStrategy;

std::string sides[3] = {"Nx", "Ny", "Nz"};

NuFFT_OMP::NuFFT_OMP () {

}

NuFFT_OMP::~NuFFT_OMP () {

	//for (int i = 0; i < NTHREADS; i++)
	//	nfft::finalize (&m_fplan[i], &m_iplan[i]);

	delete [] m_N;
	delete [] m_n;
}

RRSModule::error_code
NuFFT_OMP::Init () {

	RRSModule::error_code error = OK; 

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
	// --------------------------------------

	Attribute("maxit",   &m_maxit);
	Attribute("epsilon", &m_epsilon);
	Attribute("verbose", &m_verbose);

	int      m           = 1;
	double   alpha       = 1.0;

	Attribute("m",       &m);
	Attribute("alpha",   &alpha);

   	m_n[0] = ceil (m_N[0]*alpha);
	m_n[1] = ceil (m_N[1]*alpha);

	// Number of samples
	m_shots = m_raw.Dim(LIN);
	m_M     = m_raw.Dim(COL);

	// Initialise FT plans ------------------
	
	for (int i = 0; i < NTHREADS; i++)
		nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan[i], &m_iplan[i], m_epsilon);
	// --------------------------------------

	return error;

}

RRSModule::error_code
NuFFT_OMP::Process () {

	// Some variables
	RRSModule::error_code error = OK;

	printf ("Processing NuFFT_OMP ...\n");
	ticks start = getticks();

	// Kspace adjustment (Don't know yet why necessary)
	m_kspace = m_kspace / (1/(GAMMA/128*m_N[0]));

	int imgsize = 2;
	for (int i = 0; i < m_dim; i++)
		imgsize *= m_N[i];

	Matrix <raw> tmp;

	for (int i = 0; i < m_dim; i++)
		tmp.Dim(i) = m_N[i];

	// Store all arms separately?
	if (m_verbose)
		tmp.Dim(m_dim) = m_shots + 1;

	tmp.Reset();

	double* m_ftout    = (double*) malloc (imgsize * m_shots * sizeof(double)); 

	/*#pragma omp parallel default (shared) 
	{
		
	omp_set_num_threads(NTHREADS);*/
		int tid      = 0;//omp_get_thread_num();
		
		double* m_ftin     = (double*) malloc (2     * m_M * sizeof(double)); 
		//double* m_ftk      = (double*) malloc (m_dim * m_M * sizeof(double)); 
		//double* m_ftw      = (double*) malloc (        m_M * sizeof(double)); 

		//#pragma omp for

		for (int j = 0; j < m_shots; j++) {
			
			int     os         = j * m_M;

			// Copy data from incoming matrix to the nufft input array
			for (int i = 0; i < m_M; i++) {
				m_ftin[2*i  ] = (m_raw[i + os]).real();
				m_ftin[2*i+1] = (m_raw[i + os]).imag();
			}

			// Copy k-space and weights to allocated memory
			memcpy (&(m_fplan[tid].x), &m_kspace[os * m_dim], m_M * m_dim * sizeof(double));
			memcpy (&(m_iplan[tid].w), &m_helper[os]        , m_M *         sizeof(double));
			printf ("Made %i\n", 0);

			// Precompute PSI
			nfft::weights (&m_fplan[tid], &m_iplan[tid]);
			printf ("Made %i\n", 1);
			
			//nfft::ift     (&m_fplan[tid], &m_iplan[tid], m_ftin, &m_ftout[j * imgsize], m_maxit, m_epsilon);
			printf ("Made %i\n", 2);
			

			if (m_verbose)
				for (int i = 0; i < imgsize/2; i++)
					tmp[(j+1) * imgsize/2 + i] = raw(m_ftout[2*i + j*imgsize], m_ftout[2*i+1 + j*imgsize]);
			printf ("Made %i\n", 1);
			
			
		}

		//free (m_ftk);
		//free (m_ftw);
		free (m_ftin);

		int chunk    = imgsize/2/NTHREADS; 
		
		//#pragma omp for schedule(dynamic,chunk)
		
		for (int i = 0; i < imgsize/2; i++) {
			for (int j = 0; j < m_shots; j++)
				tmp[i] += raw(m_ftout[2*i + j*imgsize], m_ftout[2*i+1 + j*imgsize]); 
		
	}
	
	printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), start) / ClockRate());
	
	free (m_ftout);

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

