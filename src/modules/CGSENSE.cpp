#include "CGSENSE.hpp"
#include "nfftstub.h"

using namespace RRStrategy;

CGSENSE::CGSENSE () {

	int   d      = 3;
	m_N = new int[3];
	m_n = new int[3];

	FILE* fin;

	bool     weight      = true;

	unsigned infft_flags = CGNR | PRECOMPUTE_DAMP;  /* Flags for the infft */
	int      m           = 6;
	double   alpha       = 2.0;

	m_config->Attribute ("nx",      &m_N[0]);
	m_config->Attribute ("ny",      &m_N[1]);
	m_config->Attribute ("ny",      &m_N[2]);
	m_config->Attribute ("epsilon", &m_epsilon);
	m_config->Attribute ("maxit",   &m_maxit);
	m_config->Attribute ("cgconv",  &m_cgconv);
	m_config->Attribute ("ndim",    &d);
	
	nfft::init (d, m_N, m_helper.Dim(COL), m_n, m, &m_fplan, &m_iplan, m_epsilon);

	/* Copy incoming data to m_temp             */
	Matrix < raw > m_temp = (m_raw);

	// Reshape for outgoing format
	Attribute ("Nx",             &m_raw.Dim(COL));
	Attribute ("Ny",             &m_raw.Dim(LIN));
	Attribute ("Nz",             &m_raw.Dim(SLC));

	for (int i = 3; i < INVALID_DIM; i++)
		m_raw.Dim(i) = 1;

	// Keep iterations [LOTS OF DATA]? Allocate data.
	Attribute ("verbose", (int*) &m_verbose);
	if (m_verbose)
		m_raw.Dim(SET) = m_iter;

	m_raw.Reset();
	
	// Maximum k-vector reached
	kmax[0] = 0.0;
	for (int i = 0; i < m_helper.Dim(COL); i++)
		if (kmax[0] < m_helper(i,0,0))
			kmax[0] = m_helper(i,0,0);
	kmax[1] = 0.0;
	for (int j = 0; j < m_helper.Dim(LIN); j++)
		if (kmax[1] < m_helper(0,j,0))
			kmax[1] = m_helper(0,j,0);
	kmax[2] = 0.0;
	for (int k = 0; k < m_helper.Dim(SLC); k++)
		if (kmax[2] < m_helper(0,0,k))
			kmax[2] = m_helper(0,0,k);
	

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
	out->Reset();
	
	// Some dimensions
	int        ncoils   = sm->Dim(CHA);
	int        nsamples = out->Size(); 
	int        ndim     = in->Dim(COL);

	// Create container for FT input
	double     ftin  [2 *  in->Size()];

	// Loop over coils, Elementwise multiplication of maps with in (s.*in), ft and store in out
	for (int j = 0; j < ncoils; j++) {

		double ftout [2 * out->Size()/ncoils];
		int    pos = j*sm->Dim(COL)*sm->Dim(LIN);
	
		for (int i = 0; i < in->Size(); i++) {
			
			raw tmp = sm->at(pos + i) * in->at(i);
			ftin[2*i  ] = tmp.real(); 
			ftin[2*i+1] = tmp.imag(); 

		}

		nfft::ft (np, ftin, ftout);

		for (int i = 0; i < out->Size()/ncoils; i++) 
			out->at(pos + i) = raw(ftout[2*i],ftout[2*i+1]);

	}

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
	out->Reset();

	// Some dimensions
	int         ncoils   = sm->Dim(CHA);
	int         nsamples = in->Size(); 
	int         ndim     = out->Dim(COL);

	// Container for FT input
	double      ftin [2 *  in->Size()/ncoils]; 

	// Loop over coils, Inverse FT every signal in *in, 
	// Sum elementwise mutiplied images with according sensitivity maps 
	for (int j = 0; j < ncoils; j++) {
		
		double  ftout[2 * out->Size()];
		int     pos = j*in->Dim(COL)*in->Dim(LIN);

		for (int i = 0; i < in->Size()/ncoils; i++) {
			ftin[2*i  ] = (in->at(pos + i)).real();
			ftin[2*i+1] = (in->at(pos + i)).imag();
		}
		
		nfft::ift (np, spc, ftin, ftout, maxit, epsilon);

		for (int i = 0; i < out->Size(); i++) {
			raw sens = sm->at(pos + i);
			out->at(i) += raw(ftout[2*i] * sens.real(), ftout[2*i+1] * -sens.imag());
		}

	}
	
	return OK;
	
}



RRSModule::error_code
CGSENSE::Process () {
	
	Matrix<raw> a, b, p, q, r, r_new;
	
	EH (&m_temp, &m_sens, &m_fplan, &m_iplan, m_epsilon, m_maxit, &a);
	p = a;
	r = a;
	b.Dim(COL) = a.Dim(COL);
	b.Dim(LIN) = a.Dim(LIN);
	b.Reset();

	double residue = 1.0;
	double delta   = 1.0;
	
	// CG iterations
	for (int iter = 0; iter < m_iter; iter++) {
		
		float       rn;
		float       an;
		float       rnewn;
		raw         rtmp;
		Matrix<raw> mtmp;

		rn = r.norm().real();
		an = a.norm().real();

		delta = rn/an;

		if (delta < m_cgconv)
			break;
		
		// q  = eh(e(p , sensitivity, k), sensitivity, k);
		E  (&p,    &m_sens, &m_fplan,                               &mtmp);
		EH (&mtmp, &m_sens, &m_fplan, &m_iplan, m_epsilon, m_maxit, &q);

		// b  = b + r(:)'*r(:)/(p(:)'*q(:))*p;
		rtmp  = (rn / (p.dotc(q)));
		mtmp  = p * rtmp;
		b     = b + mtmp;
		
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

