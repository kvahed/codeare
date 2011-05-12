#include "CGSENSE.hpp"
#include "nfftstub.h"

#include <vector>

using namespace RRStrategy;

CGSENSE::CGSENSE () {

}


RRSModule::error_code 
CGSENSE::Init() {
	
	RRSModule::error_code error = OK; 

	// Verbosity ----------------------------

	Attribute ("verbose", &m_verbose);
	// --------------------------------------

	// Dimensions ---------------------------
	Attribute("dim",     &m_dim);

	m_N = new int[m_dim];
	m_n = new int[m_dim];

	Attribute ("Nx",      &m_N[0]);
	Attribute ("Ny",      &m_N[1]);
	Attribute ("M",       &m_M);

	m_M = m_helper.Size();
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
	int      m           = 6;
	double   alpha       = 2.0;

   	m_n[0] = ceil (m_N[0]*alpha);
	m_n[1] = ceil (m_N[1]*alpha);
	// --------------------------------------

	// Initialise FT plans ------------------
	nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan, &m_iplan, m_epsilon);
	// --------------------------------------

	// Allocate RAM for fix size memories (weights and k-space)
	m_ftk      = (double*) malloc (    m_dim  * m_M    * sizeof(double)); 
	m_ftw      = (double*) malloc (             m_M    * sizeof(double)); 
	// --------------------------------------

	return error;

}


/**
 * @brief               Compute left hand side (i.e. Multiply E with spatial (image) data)
 *                      Forward NFFT in and elementwise multiply with spatial sensitivity of every channel 
 * 
 * @param  in           Original discretised sample O (Nx x Ny)
 * @param  sm           Sensitivity maps            O (Nx x Ny x Nc)
 * @param  np           Non-Cartesian strategy for non uniform ft
 * @param  out          Result                      O (Nk x Nc)
 */
RRSModule::error_code 
E  (Matrix<raw>* in, Matrix<raw>* sm, nfft_plan* np, Matrix<raw>* out) {

	// Clear output container
	out->Zero();
	
	// Some dimensions
	int        ncoils   = sm->Dim(CHA);
	int        nsamples = out->Size() / ncoils;
	int        imgsize  = in->Size();

	// Create container for FT input
	double*    ftin     = (double*) malloc (2 * imgsize  * sizeof(double));
	double*    ftout    = (double*) malloc (2 * nsamples * sizeof(double));

	// Loop over coils, Elementwise multiplication of maps with in (s.*in), ft and store in out
	for (int j = 0; j < ncoils; j++) {

		int     ipos     = j * imgsize;
		int     spos     = j * nsamples;
	
		// Copy data to FT
		for (int i = 0; i < imgsize; i++) {
			
			raw tmp     = sm->at(ipos + i) * in->at(i);

			ftin[2*i  ] = tmp.real(); 
			ftin[2*i+1] = tmp.imag(); 

		}

		// Forward ft
		nfft::ft (np, ftin, ftout);

		// Copy FTed data back
		for (int i = 0; i < nsamples; i++) 
			out->at(spos + i) = raw(ftout[2*i],ftout[2*i+1]);

	}

	// Free RAM
	free (ftout);
	free (ftin);

	// Return success
	return OK;

}


/**
 * @brief               Compute right hand side (i.e. Multiply E^H, Hermitian counterpart to E, with k-space data)
 *
 * @param  in           K-space samples along trajectory O (Nk x Nc)
 * @param  sm           Sensitivity maps                 O (Nx x Ny x Nc)
 * @param  np           NuFFT plan
 * @param  spc          Solver plan
 * @param  epsilon      Convergence criterium for ift (default 3e-7)
 * @param  maxit        Maximum number of solver iterations (default 3)
 * @param  out          Returned product                 O (Nx x Ny)
 */
RRSModule::error_code
EH (Matrix<raw>* in, Matrix<raw>* sm, nfft_plan* np, solver_plan_complex* spc, double epsilon, int maxit, Matrix<raw>* out) {

	// Clear outgoing container
	out->Zero();

	// Some dimensions
	int        ncoils   = sm->Dim(CHA);
	int        nsamples = in->Size() / ncoils;
	int        imgsize  = out->Size();

	// Containers for FT I/O
	double*     ftin     = (double*) malloc (2 * nsamples * sizeof(double));
	double*     ftout    = (double*) malloc (2 * imgsize  * sizeof(double));

	// Loop over coils, Inverse FT every signal in *in, 
	// Sum elementwise mutiplied images with according sensitivity maps 
	for (int j = 0; j < ncoils; j++) {
		int    spos      = j * nsamples;
		int    ipos      = j * imgsize;

		// Copy to iFT
		for (int i = 0; i < nsamples; i++) {
			ftin[2*i  ] = (in->at(spos + i)).real();
			ftin[2*i+1] = (in->at(spos + i)).imag();
		}
		
		// Inverse FT
		nfft::ift (np, spc, ftin, ftout, maxit, epsilon);

		// Copy back from iFT
		for (int i = 0; i < out->Size(); i++) {
			raw sens = sm->at(ipos + i);
			out->at(i) += raw(ftout[2*i], ftout[2*i+1]) * conj(sens);
		}

	}
	
	// Free RAM
	free (ftout);
	free (ftin);

	return OK;
	
}


RRSModule::error_code
CGSENSE::Process () {
	
	RRSModule::error_code error = CG_DID_NOT_CONVERGE;
	static clock_t runtime = clock();

	Matrix<raw> a, p, q, r, r_new;

	for (int i = 0; i < INVALID_DIM; i++)
		a.Dim(i) = 1;
	for (int i = 0; i < m_dim      ; i++)
		a.Dim(i) = m_N[i];
	a.Reset();

	m_raw = m_raw * 10000;

	m_kspace = m_kspace / (1/(GAMMA/128*m_N[0]));

	double*     ftk      = (double*) malloc (   m_kspace.Size() * sizeof(double)); 
	double*     ftw      = (double*) malloc (   m_helper.Size() * sizeof(double)); 

	memcpy (ftw, &m_helper[0], m_helper.Size()*sizeof(double));
	memcpy (ftk, &m_kspace[0], m_kspace.Size()*sizeof(double));
	
	nfft::kspace  (&m_fplan,           ftk);
	nfft::weights (&m_fplan, &m_iplan, ftw);

	EH (&m_raw, &m_rhelper, &m_fplan, &m_iplan, m_epsilon, m_maxit, &a);
	p = Matrix<raw>(a);
	r = Matrix<raw>(a);

	// Prepare q
	q.Dim(COL) = a.Dim(COL);
	q.Dim(LIN) = a.Dim(LIN);
	q.Reset();

	// Out going image
	// Resize m_raw for output
	for (int i = 0; i < INVALID_DIM; i++)
		m_raw.Dim(i) = 1;
	for (int i = 0; i < m_dim; i++)
		m_raw.Dim(i) = m_N[i];
	m_raw.Reset();

	std::vector<double> res;
	
	Matrix<raw> sigtmp;
	sigtmp.Dim (COL) = m_M;
	sigtmp.Dim (LIN) = 8;
	sigtmp.Reset();

	Matrix<raw> imgtmp;
	imgtmp.Dim (COL) = m_N[0];
	imgtmp.Dim (LIN) = m_N[1];
	imgtmp.Reset();

	float       rn    = 0.0;
	float       an    = 0.0;
	float       rnewn = 0.0;
	raw         rtmp  = raw(0.0,0.0);
	
	// CG iterations
	for (int i = 0; i < m_cgmaxit; i++) {

		rn = r.norm().real();
		an = a.norm().real();
		
		res.push_back(rn/an);
		
		printf ("%i: CG residuum: %.9f\n", i, res.at(i));
		
		if (res.at(i) < m_cgeps) {
			error = OK;
			break;
		}
		
		// q  = eh(e(p , sensitivity, k), sensitivity, k);

		E  (&p,      &m_rhelper, &m_fplan,                               &sigtmp);
		EH (&sigtmp, &m_rhelper, &m_fplan, &m_iplan, m_epsilon, m_maxit, &q);
		
		// b     = b + r(:)'*r(:)/(p(:)'*q(:))*p;
		rtmp     = (rn / (p.dotc(q)));
		imgtmp   = p * rtmp;
		m_raw    = m_raw + imgtmp;
		
		// r_new = r - r(:)'*r(:)/(p(:)'*q(:))*q; 
		imgtmp   = q * rtmp;
		r_new    = r - imgtmp;
		
		// p     = r_new + r_new(:)'*r_new(:)/(r(:)'*r(:))*p
		rnewn    = r_new.norm().real();
		rtmp     = rnewn/rn;
		imgtmp   = p * rtmp;
		p        = r_new + imgtmp;
		
		// r     = r_new
		r        = r_new;

	}

	if (m_verbose)
		a.dump("share/cgsense/test.h5");

	free (ftk);
	free (ftw);

	runtime = clock() - runtime;
	printf ("Processing NuFFT took: %.4f seconds.\n", runtime / 1000000.0);

	return error;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

