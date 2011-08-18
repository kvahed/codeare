#include "CGSENSE.hpp"
#include "nfftstub.h"
#include "Noise.hpp"
#include "SEM.hpp"
#include <math.h>

#include <vector>

std::string sides[3] = {"Nx", "Ny", "Nz"};

using namespace RRStrategy;


CGSENSE::~CGSENSE () {
	this->Finalise();
}


RRSModule::error_code 
CGSENSE::Finalise () {

	if (m_initialised)
		for (int i = 0; i < NTHREADS || i < m_Nc; i++)
			nfft::finalize (&m_fplan[i], &m_iplan[i]);

	delete [] m_N;
	delete [] m_n;
	return OK;

}


RRSModule::error_code 
CGSENSE::Init() {

	printf ("Intialising CG-SENSE ...\n");

	RRSModule::error_code error = OK; 

	m_initialised = false;

	// Image space dimensions ----------------
	m_N   = new int[3];

	// Oversampling --------------------------
	m_n   = new int[3];

	// Some defaults
	for (int i = 0; i < 3; i++) {
		m_N[i] = 1; 
		m_n[i] = 0;
	}

	m_testcase = 0;
	m_verbose  = 0;
	m_noise    = 0;
	m_dim      = 1;
	m_M        = 0;

	// Dimensions ---------------------------

	Attribute("dim",       &m_dim);
	printf ("  dimensions: %iD \n", m_dim);


	if (m_dim < 2 || m_dim > 3) {
		printf ("%s only supports 2-3 dimensions. %i was specified.\n", Name(), m_dim);
		return UNSUPPORTED_DIMENSION;
	}

	for (int i = 0; i < m_dim; i++)
		Attribute (sides[i].c_str(),       &m_N[i]);

	for (int i = 0; i < m_dim; i++)
		if (m_N[i] < 16) {
			printf ("%s only supports image matrix sides >= 16. (%ix%ix%i) was specified.\n", Name(), m_N[0], m_N[1], m_N[2]);
			return UNSUPPORTED_IMAGE_MATRIX;
		}

	Attribute ("M",         &m_M);

	if (m_M == 0) {
		printf ("Initialising %s with %i mesurement nodes? Check configuration! FAILED: Bailing out!\n", Name(), m_M);
		return ZERO_NODES;
	}

	Attribute ("Nc",        &m_Nc);

	if (m_Nc == 0) {
		printf ("Initialising %s with %i channels? Check configuration! FAILED: Bailing out!\n", Name(), m_M);
		return CGSENSE_ZERO_CHANNELS;
	}

	// --------------------------------------

	// Verbosity ----------------------------

	Attribute ("testcase",  &m_testcase);
	// --------------------------------------

	// Verbosity ----------------------------

	Attribute ("verbose",   &m_verbose);
	printf ("  verbose feedback: %i \n", m_verbose);
	// --------------------------------------

	// Noise --------------------------------

	Attribute ("noise",   &m_noise);
	printf ("  gaussian white noise (normalised): %.9f \n", m_noise);
	// --------------------------------------

	// CG convergence and break criteria ----

	Attribute ("cgeps",   &m_cgeps);
	Attribute ("cgmaxit", &m_cgmaxit);
	printf ("  maximum #iterations: %i \n", m_cgmaxit);
	printf ("  convergence criterium: %.9f \n", m_cgeps);
	// --------------------------------------

	// iNFFT convergence and break criteria -

	Attribute ("maxit",   &m_maxit);
	Attribute ("epsilon", &m_epsilon);
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
	
	printf ("  intialising nfft::init (%i, {%i, %i, %i}, %i, {%i, %i, %i}, %i, *, *, %.9f)\n", 
			m_dim, 
			m_N[0], m_N[1], m_N[2],
			m_M,
			m_n[0], m_n[1], m_n[2],
			m,
			m_epsilon);

	for (int i = 0; i < NTHREADS || i < m_Nc; i++)
		nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan[i], &m_iplan[i], m_epsilon);
	// --------------------------------------

	m_initialised = true;

	printf ("... done.\n\n");

	return error;

}


RRSModule::error_code
CGSENSE::Process () {

	RRSModule::error_code error = OK;
	
	// CG matrices ----------------------------------------------------
	Matrix<raw> p, q, r;

	// Add white noise? (Only for testing) ----------------------------
	if (m_noise > 0.0)
		AddPseudoRandomNoise (&m_raw, (float)m_noise);
	
	// ----------------------------------------------------------------
	for (int i = 0; i < m_dim      ; i++)
		p.Dim(i) = m_N[i];
	p.Reset();
	
	// Set k-space and weights ----------------------------------------
	for (int i = 0; i < NTHREADS || i < m_Nc; i++) {
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

	// Create test data (Incoming data is image space) ----------------
	if (m_testcase) {
		E  (&m_rhelper, &m_sens, m_fplan, &stmp, m_dim);
		m_rhelper = stmp;
		m_raw     = stmp;
	}


	// Start CG routine and runtime -----------------------------------
	ticks cgstart = getticks();

	// First left side action -----------------------------------------
	EH (&m_raw, &m_sens, m_fplan, m_iplan, m_epsilon, m_maxit, &p, m_dim);

	r = p;
	q = p;

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

	an = pow(p.Norm().real(),2.0);

	m_rhelper.Zero();

	for (int i = 0; i < m_cgmaxit; i++, iters++) {

		rn = pow(r.Norm().real(),2.0);

		res.push_back(rn/an);

		printf ("  %03i: CG residuum: %.9f\n", i, res.at(i));

		// Convergence ? ----------------------------------------------
		if (std::isnan(res.at(i)) || res.at(i) <= m_cgeps)
			break;
		
		// CG step ----------------------------------------------------
		E  (&p,    &m_sens, m_fplan,                              &stmp, m_dim);
		EH (&stmp, &m_sens, m_fplan, m_iplan, m_epsilon, m_maxit, &q   , m_dim);

		rtmp      = (rn / (p.dotc(q)));
		m_raw     = (p * rtmp) + m_raw;
		stmp     *= rtmp;
		m_rhelper = m_rhelper + stmp;
		r         = - (q * rtmp) + r ;
		rnewn     = pow(r.Norm().real(),2.0);
		rtmp      = rnewn/rn;
		p         = (p * rtmp) + r ;

		// Verbose out put keeps all intermediate steps ---------------
		if (m_verbose) {
			memcpy (&istore[i *     m_raw.Size()],     &m_raw[0],     m_raw.Size() * sizeof(double));
			memcpy (&sstore[i * m_rhelper.Size()], &m_rhelper[0], m_rhelper.Size() * sizeof(double));
		}

		if (std::isnan(res.at(i)))
			break;

	}

	// Report timimng -------------------------------------------------
	printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());

	// Verbose output needs to 
	if (m_verbose) {

		// All intermediate images ------
		m_raw.Dim(m_dim) = iters;
		m_raw.Reset();
		memcpy (    &m_raw[0], &istore[0],     m_raw.Size() * sizeof(double));

		// Pulses (Excitation) ----------
		m_rhelper.Dim(CHA) = iters;
		m_rhelper.Reset();
		memcpy (&m_rhelper[0], &sstore[0], m_rhelper.Size() * sizeof(double));

		// CG residuals ------------------
		m_helper.Dim(COL) = iters;
		for (int i = 1; i < INVALID_DIM; i++)
			m_helper.Dim(i) = 1;
		m_helper.Reset();

		for (int i = 0; i < iters; i++)
			m_helper[i] = res.at(i);

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


