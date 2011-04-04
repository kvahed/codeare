#include "CGSENSE.h"

using namespace RRStrategy;

CGSENSE::CGSENSE () {

	m_N = new int[2];
	m_n = new int[2];

	FILE* fin;

	bool     weight      = true;

	unsigned infft_flags = CGNR | PRECOMPUTE_DAMP;  /* Flags for the infft */
	int      m           = 6;
	double   alpha       = 2.0;

	m_config->Attribute ("nx",      &m_N[0]);
	m_config->Attribute ("ny",      &m_N[1]);
	m_config->Attribute ("epsilon", &m_epsilon);
	m_config->Attribute ("maxit",   &m_maxit);
	m_config->Attribute ("cgconv",  &m_cgconv);
	

	/* initialise my_plan */
	m_n[0] = ceil(m_N[0]*alpha);
	m_n[1] = ceil(m_N[1]*alpha);

	nfft_init_guru (&m_plan, 2, m_N, m_raw.Dim(COL), m_n, m, 
					PRE_PHI_HUT| PRE_PSI| MALLOC_X| MALLOC_F_HAT| MALLOC_F|	FFTW_INIT| FFT_OUT_OF_PLACE, FFTW_MEASURE| FFTW_DESTROY_INPUT);

	/* Precompute lin psi if set */
	if(m_plan.nfft_flags & PRE_LIN_PSI)
		nfft_precompute_lin_psi(&m_plan);

	/* Set the flags for the infft*/
	if (weight)
		infft_flags = infft_flags | PRECOMPUTE_WEIGHT;

	/* Get the weights */
	if (m_iplan.flags & PRECOMPUTE_WEIGHT)
		for(int j = 0; j < m_plan.M_total; j++)
			memcpy(&m_iplan.w[0], &m_weights[0], sizeof(double)*m_weights.Size());

	/* initialise my_iplan, advanced */
	solver_init_advanced_complex (&m_iplan, (mv_plan_complex*)&m_plan, infft_flags);

	/* get the damping factors */
	if (m_iplan.flags & PRECOMPUTE_DAMP)
		for (int j = 0; j < m_N[0]; j++) {
			for (int k = 0; k < m_N[1]; k++) {
				int    j2 = j - m_N[0]/2;
				int    k2 = k - m_N[1]/2;
				double r  = sqrt (j2*j2 + k2*k2);
				if (r > (double)m_N[0]/2)
					m_iplan.w_hat[j*m_N[0]+k] = 0.0;
				else
					m_iplan.w_hat[j*m_N[0]+k] = 1.0;
			}
		}

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
 *
 * @param  in           Original discretised sample O (Nx x Ny)
 * @param  sm           Sensitivity maps            O (Nx x Ny x Nc)
 * @param  ncs          Non-Cartesian strategy for non uniform ft
 * @param  out          Result                      O (Nk x Nc)
 */
RRSModule::error_code 
E  (Matrix<raw>* in, Matrix<raw>* sm, nfft_plan* plan, Matrix<raw>* out) {

	int    ncoils   = sm->Dim(CHA);
	int    nsamples = out->Size(); 
	int    ndim     = in->Dim(COL);
	
	for (int j = 0; j < ncoils; j++) {

		
		
	}


	return OK;

}


/**
 * @brief               Compute right hand side (i.e. Multiply E^H, Hermitian counterpart to E, with k-space data)
 *
 * @param  in           K-space samples along trajectory O (Nk x Nc)
 * @param  sm           Sensitivity maps                 O (Nx x Ny x Nc)
 * @param  ncs          Non-Cartesian strategy for non uniform ft
 * @param  out          Returned product                 O (Nx x Ny)
 */
RRSModule::error_code
EH (Matrix<raw>* in, Matrix<raw>* sm, solver_plan_complex* plan, double epsilon, int maxit, Matrix<raw>* out) {

	int    ncoils   = sm->Dim(CHA);
	int    nsamples = in->Size(); 
	int    ndim     = out->Dim(COL);

	for (int j = 0; j < ncoils; j++) {
		
		for (int i=0; i < nsamples; i++) {
			plan->y[i][0] = (in->channel(j))[i].real();
			plan->y[i][1] = (in->channel(j))[i].imag();
		}
		
		double t        = nfft_second();

		/* inverse trafo */
		solver_before_loop_complex(plan);
		
		for (int i = 0; i < maxit; i++)  {

			/* break if dot_r_iter is smaller than epsilon*/
			if(plan->dot_r_iter<epsilon)
				break;

#ifdef VERBOSE
			fprintf(stderr,"%e,  %i of %i\n",sqrt(plan->dot_r_iter), i+1, maxit);
#endif

			solver_loop_one_step_complex(plan);

		}
		
		t = nfft_second()-t;
		
		for (int i = 0; i < out->Size(); i++)
			out->at(i) += sm->at(i) * raw(plan->f_hat_iter[i][0], plan->f_hat_iter[i][1]);
		
	}
	
	return OK;
	
}



RRSModule::error_code
CGSENSE::Process () {
	
	Matrix<raw> a, b, p, q, r, r_new;
	
	EH (&m_temp, &m_sens, &m_iplan, m_epsilon, m_maxit, &a);
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
		E  (&p,    &m_sens, &m_fplan,                     &mtmp);
		EH (&mtmp, &m_sens, &m_iplan, m_epsilon, m_maxit, &q);

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

