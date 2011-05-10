#include "CGSENSE.hpp"
#include "nfftstub.h"

#include <vector>

using namespace RRStrategy;

CGSENSE::CGSENSE () {

}


RRSModule::error_code 
CGSENSE::Init() {
	
	RRSModule::error_code error = OK; 

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

		int     pos     = j * imgsize;
	
		// Copy data to FT
		for (int i = 0; i < imgsize; i++) {
			
			raw tmp     = sm->at(pos + i) * in->at(i);

			ftin[2*i  ] = tmp.real(); 
			ftin[2*i+1] = tmp.imag(); 

		}

		// Forward ft
		nfft::ft (np, ftin, ftout);

		// Copy FTed data back
		for (int i = 0; i < nsamples; i++) 
			out->at(pos + i) = raw(ftout[2*i],ftout[2*i+1]);

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
		
		int     pos      = j * nsamples;

		// Copy to iFT
		for (int i = 0; i < nsamples; i++) {
			ftin[2*i  ] = (in->at(pos + i)).real();
			ftin[2*i+1] = (in->at(pos + i)).imag();
		}
		
		// Inverse FT
		nfft::ift (np, spc, ftin, ftout, maxit, epsilon);

		// Copy back from iFT
		for (int i = 0; i < out->Size(); i++) {
			raw sens = sm->at(pos + i);
			out->at(i) += raw(ftout[2*i] * sens.real(), ftout[2*i+1] * -sens.imag());
		}

	}
	
	// Free RAM
	free (ftout);
	free (ftin);

	return OK;
	
}


RRSModule::error_code
CGSENSE::Process () {
	
	Matrix<raw> a, p, q, r, r_new;

	// First application of the right hand side
	EH (&m_raw, &m_sens, &m_fplan, &m_iplan, m_epsilon, m_maxit, &a);
	p = a;
	r = a;

	// Out going image
	// Resize m_raw for output
	for (int i = 0; i < INVALID_DIM; i++)
		m_raw.Dim(i) = 1;
	for (int i = 0; i < m_dim; i++)
		m_raw.Dim(i) = m_N[i];
	m_raw.Reset();

	std::vector<double> res;
	
	// CG iterations
	for (int i = 0; i < m_cgmaxit; i++) {
		
		float       rn;
		float       an;
		float       rnewn;
		raw         rtmp;
		Matrix<raw> mtmp;
		
		rn = r.norm().real();
		an = a.norm().real();
		
		res.push_back(rn/an);
		
		if (res.at(i) < m_cgeps)
			break;
		
		// q  = eh(e(p , sensitivity, k), sensitivity, k);
		E  (&p,    &m_sens, &m_fplan,                               &mtmp);
		EH (&mtmp, &m_sens, &m_fplan, &m_iplan, m_epsilon, m_maxit, &q);
		
		// b  = b + r(:)'*r(:)/(p(:)'*q(:))*p;
		rtmp  = (rn / (p.dotc(q)));
		mtmp  = p * rtmp;
		m_raw = m_raw + mtmp;
		
		// r_new = r - r(:)'*r(:)/(p(:)'*q(:))*q; 
		mtmp  = q * rtmp;
		r_new = r - mtmp;
		
		// p  = r_new + r_new(:)'*r_new(:)/(r(:)'*r(:))*p
		rnewn = r_new.norm().real();
		rtmp  = rnewn/rn;
		mtmp  = p * rtmp;
		p     = r_new + mtmp;
		
		// r  = r_new
		r     = r_new;
		
	}
	
	return OK;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

