#include "CGSENSE.hpp"
#include "nfftstub.h"
#include "Noise.hpp"
#include "SEM.hpp"

#include <vector>

std::string sides[3] = {"Nx", "Ny", "Nz"};

using namespace RRStrategy;


CGSENSE::~CGSENSE () {

	free (m_ftw);
	free (m_ftk);

	for (int i = 0; i < NTHREADS; i++)
		nfft::finalize (&m_fplan[i], &m_iplan[i]);

	delete [] m_N;
	delete [] m_n;

};


RRSModule::error_code 
CGSENSE::Init() {

	RRSModule::error_code error = OK; 

	m_N   = new int[3];
	m_n   = new int[3];

	for (int i = 0; i < 3; i++) {
		m_N[i] = 0; 
		m_n[i] = 0;
	}

	m_testcase = 0;
	m_verbose  = 0;
	m_noise    = 0;

	// Verbosity ----------------------------

	Attribute ("testcase",  &m_testcase);
	// --------------------------------------

	// Verbosity ----------------------------

	Attribute ("verbose",   &m_verbose);
	// --------------------------------------

	// Noise --------------------------------

	Attribute ("noise",   &m_noise);
	// --------------------------------------

	// Dimensions ---------------------------

	Attribute("dim",       &m_dim);

	for (int i = 0; i < m_dim; i++)
		Attribute (sides[i].c_str(),       &m_N[i]);

	m_M   = m_helper.Size();
	// --------------------------------------

	// iNFFT convergence and break criteria -

	Attribute ("maxit",   &m_maxit);
	Attribute ("epsilon", &m_epsilon);
	// --------------------------------------

	// CG convergence and break criteria

	Attribute ("cgeps",   &m_cgeps);
	Attribute ("cgmaxit", &m_cgmaxit);
	// --------------------------------------

	// Oversampling -------------------------

	int      m           = 0;
	double   alpha       = 0.0;

	Attribute ("m",       &m);
	Attribute ("alpha",   &alpha);

	for (int i = 0; i < m_dim; i++)
		m_n[i] = ceil (m_N[i]*alpha);
	// --------------------------------------

	// Initialise FT plans ------------------
	
	for (int i = 0; i < NTHREADS; i++)
		nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan[i], &m_iplan[i], m_epsilon);
	// --------------------------------------

	// Allocate RAM for fix size memories (weights and k-space)

	m_ftk      = (double*) malloc (    m_dim  * m_M    * sizeof(double)); 
	m_ftw      = (double*) malloc (             m_M    * sizeof(double)); 
	// --------------------------------------

	return error;

}


RRSModule::error_code
CGSENSE::Process () {

	RRSModule::error_code error = OK;

	Matrix<raw> a, p, q, r, r_new;

	if (m_noise > 0.0)
		AddPseudoRandomNoise (&m_raw, (float)m_noise);
	
	for (int i = 0; i < INVALID_DIM; i++)
		a.Dim(i) = 1;
	for (int i = 0; i < m_dim      ; i++)
		a.Dim(i) = m_N[i];

	a.Reset();
	
	memcpy (m_ftw, &m_helper[0], m_helper.Size()*sizeof(double));
	memcpy (m_ftk, &m_kspace[0], m_kspace.Size()*sizeof(double));
	
	for (int i = 0; i < NTHREADS; i++) {
		nfft::kspace  (&m_fplan[i],              m_ftk);
		nfft::weights (&m_fplan[i], &m_iplan[i], m_ftw);
	}

	m_sens = m_rhelper;

	Matrix<raw> store;
	for (int i = 0; i < INVALID_DIM; i++)
		store.Dim(i) = 1;

	if (m_verbose == 1) 
		for (int i = 0; i < m_dim; i++)
			store.Dim(i) = m_N[i];

	store.Dim(m_dim) = m_cgmaxit;
	store.Reset();

	m_rhelper = m_raw;

	Matrix<raw> sigtmp;
	sigtmp.Dim (COL) = m_M;
	sigtmp.Dim (LIN) = 8;
	sigtmp.Reset();

	Matrix<raw> imgtmp;
	
	for (int i = 0; i < m_dim; i++)
		imgtmp.Dim (i) = m_N[i];

	imgtmp.Reset();

	if (m_testcase) {
		E  (&m_rhelper, &m_sens, &m_fplan[0], &sigtmp, m_dim);
		m_rhelper = sigtmp;
		m_raw     = sigtmp;
	}

	ticks cgstart = getticks();

	EH (&m_raw, &m_sens, &m_fplan[0], &m_iplan[0], m_epsilon, m_maxit, &a, m_dim);
	p = a;
	r = a;
	q = a;

	// Out going image
	// Resize m_raw for output
	for (int i = 0; i < INVALID_DIM; i++)
		m_raw.Dim(i) = 1;

	for (int i = 0; i < m_dim; i++)
		m_raw.Dim(i) = m_N[i];

	m_raw.Reset();

	std::vector<double> res;
	
	float       rn    = 0.0;
	float       an    = 0.0;
	float       rnewn = 0.0;
	raw         rtmp  = raw(0.0,0.0);
	int         iters = 0;

	printf ("Processing CG-SENSE ...\n");

	// CG iterations (Pruessmann et al. (2001). MRM, 46(4), 638-51.)
	for (int i = 0; i < m_cgmaxit; i++, iters++) {

		rn = r.norm().real();
		an = a.norm().real();
		
		res.push_back(rn/an);
		
		printf ("%03i: CG residuum: %.9f\n", i, res.at(i));
		if (res.at(i) <= m_cgeps)
			break;
		
		E  (&p,      &m_sens, &m_fplan[0],                               &sigtmp, m_dim);
		EH (&sigtmp, &m_sens, &m_fplan[0], &m_iplan[0], m_epsilon, m_maxit, &q     , m_dim);
		
		rtmp      = (rn / (p.dotc(q)));
		imgtmp    = p * rtmp;
		m_raw     = m_raw + imgtmp;
		sigtmp    = sigtmp * rtmp;
		m_rhelper = m_rhelper + sigtmp;
		imgtmp    = q * rtmp;
		r_new     = r - imgtmp;
		rnewn     = r_new.norm().real();
		rtmp      = rnewn/rn;
		imgtmp    = p * rtmp;
		p         = r_new + imgtmp;
		r         = r_new;

		if (m_verbose == 1)
			memcpy (&store[i*m_raw.Size()], &m_raw[0], m_raw.Size() * sizeof(double));

	}

	printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), cgstart) / ClockRate());

	if (m_verbose == 1) {
		m_raw.Dim(2) = iters;
		m_raw.Reset();
		memcpy (&m_raw[0], &store[0], m_raw.Size() * sizeof(double));
	}

	return error;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}


