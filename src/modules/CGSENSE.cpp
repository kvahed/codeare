#include "CGSENSE.hpp"
#include "nfftstub.h"
#include "Noise.hpp"
#include "SEM.hpp"

#include <vector>

std::string sides[3] = {"Nx", "Ny", "Nz"};

using namespace RRStrategy;


CGSENSE::~CGSENSE () {

	for (int i = 0; i < NTHREADS; i++)
		nfft::finalize (&m_fplan[i], &m_iplan[i]);

	delete [] m_N;
	delete [] m_n;

}


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
	m_dim      = 1;

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

	// CG convergence and break criteria ----

	Attribute ("cgeps",   &m_cgeps);
	Attribute ("cgmaxit", &m_cgmaxit);
	// --------------------------------------

	// Oversampling -------------------------

	int      m           = 1;
	double   alpha       = 1.0;

	Attribute ("m",       &m);
	Attribute ("alpha",   &alpha);

	for (int i = 0; i < m_dim; i++)
		m_n[i] = ceil (m_N[i]*alpha);
	// --------------------------------------

	// Initialise FT plans ------------------
	
	for (int i = 0; i < NTHREADS; i++)
		nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan[i], &m_iplan[i], m_epsilon);
	// --------------------------------------

	return error;

}


RRSModule::error_code
CGSENSE::Process () {

	RRSModule::error_code error = OK;

	// CG matrices ----------------------------------------------------
	Matrix<raw> a, p, q, r, r_new;

	// Add white noise? (Only for testing) ----------------------------
	if (m_noise > 0.0)
		AddPseudoRandomNoise (&m_raw, (float)m_noise);
	
	// ----------------------------------------------------------------
	for (int i = 0; i < m_dim      ; i++)
		a.Dim(i) = m_N[i];
	a.Reset();
	
	// Set k-space and weights ----------------------------------------
	for (int i = 0; i < NTHREADS; i++) {
		memcpy (&(m_fplan[i].x[0]), &m_kspace[0], m_fplan[i].d * m_fplan[i].M_total * sizeof(double));
		memcpy (&(m_iplan[i].w[0]), &m_helper[0],                m_fplan[i].M_total * sizeof(double)) ;
		nfft::weights (&m_fplan[i], &m_iplan[i]);
	}

	// Copying sensitivities. Will use helper for Pulses --------------
	m_sens    = m_rhelper;
	m_rhelper = m_raw;

	// Out going images -----------------------------------------------
	Matrix<raw> istore;

	if (m_verbose == 1) 
		for (int i = 0; i < m_dim; i++)
			istore.Dim(i) = m_N[i];

	istore.Dim(m_dim) = m_cgmaxit;
	istore.Reset();

	// Temporary signal repository ------------------------------------ 
	Matrix<raw> stmp;
	stmp.Dim (COL) = m_M;
	stmp.Dim (LIN) = 8;
	stmp.Reset();

	// Out going signals ----------------------------------------------
	Matrix<raw> sstore = stmp;
	if (m_verbose == 1) 
		sstore.Dim (CHA) = m_cgmaxit;
	sstore.Reset();

	// Temporary imag repository --------------------------------------
	Matrix<raw> itmp;
	for (int i = 0; i < m_dim; i++)
		itmp.Dim (i) = m_N[i];
	itmp.Reset();
	
	// Create test data (Incoming data is image space) ----------------
	if (m_testcase) {
		E  (&m_rhelper, &m_sens, m_fplan, &stmp, m_dim);
		m_rhelper = stmp;
		m_raw     = stmp;
	}

	// Start CG routine and runtime -----------------------------------
	ticks cgstart = getticks();

	// First left side action -----------------------------------------
	EH (&m_raw, &m_sens, m_fplan, m_iplan, m_epsilon, m_maxit, &a, m_dim);
	p = a;
	r = a;
	q = a;

	// Out going image ------------------------------------------------
	// Resize m_raw for output
	for (int i = 0; i < INVALID_DIM; i++)
		m_raw.Dim(i) = 1;
	for (int i = 0; i < m_dim; i++)
		m_raw.Dim(i) = m_N[i];

	m_raw.Reset();
	
	// CG residuals storage and helper variables ----------------------
	std::vector<double> res;

	float       rn    = 0.0;
	float       an    = 0.0;
	float       rnewn = 0.0;
	raw         rtmp  = raw(0.0,0.0);
	int         iters = 0;

	printf ("Processing CG-SENSE ...\n");

	// CG iterations (Pruessmann et al. (2001). MRM, 46(4), 638-51.) --
	for (int i = 0; i < m_cgmaxit; i++, iters++) {

		rn = r.norm().real();
		an = a.norm().real();

		res.push_back(rn/an);
		
		printf ("%03i: CG residuum: %.9f\n", i, res.at(i));

		// Convergence ? ----------------------------------------------
		if (res.at(i) <= m_cgeps)
			break;
		
		// CG step ----------------------------------------------------
		E  (&p,    &m_sens, m_fplan,                                  &stmp, m_dim);
		EH (&stmp, &m_sens, m_fplan, m_iplan, m_epsilon, m_maxit, &q   , m_dim);
		
		rtmp      = (rn / (p.dotc(q)));
		itmp      = p * rtmp;
		m_raw     = m_raw + itmp;
		stmp      = stmp * rtmp;
		m_rhelper = m_rhelper + stmp;
		itmp      = q * rtmp;
		r_new     = r - itmp;
		rnewn     = r_new.norm().real();
		rtmp      = rnewn/rn;
		itmp      = p * rtmp;
		p         = r_new + itmp;
		r         = r_new;

		// Verbose out put keeps all intermediate steps ---------------
		if (m_verbose) {
			memcpy (&istore[i *     m_raw.Size()],     &m_raw[0],     m_raw.Size() * sizeof(double));
			memcpy (&sstore[i * m_rhelper.Size()], &m_rhelper[0], m_rhelper.Size() * sizeof(double));
		}

	}

	// Report timimng -------------------------------------------------
	printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), cgstart) / ClockRate());

	// Verbose output needs to 
	if (m_verbose) {
		m_raw.Dim(m_dim) = iters;
		m_raw.Reset();
		memcpy (    &m_raw[0], &istore[0],     m_raw.Size() * sizeof(double));

		m_rhelper.Dim(CHA) = iters;
		m_rhelper.Reset();
		memcpy (&m_rhelper[0], &sstore[0], m_rhelper.Size() * sizeof(double));
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


