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
	if (m_iplan.flags & PRECOMPUTE_WEIGHT) {
		fin = fopen ("weights.dat","r");
		for(int j = 0; j < m_plan.M_total; j++)
			fscanf (fin, "%le ", &m_iplan.w[j]);
		fclose (fin);
	}

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
 * @param  kt           k-space trajectory          O (Nk x 1)
 * @param  ncs          Non-Cartesian strategy for non uniform ft
 * @param  out          Result                      O (Nk x Nc)
 */
RRSModule::error_code 
E  (Matrix<raw>* in, Matrix<raw>* sm, Matrix<raw>* kt, nfft_plan* plan, Matrix<raw>* out) {
	
	int nc   = sm->Dim(CHA);
	int nk   = in->Dim(COL); 

	int nx   = sm->Dim(COL);
	int ny   = sm->Dim(LIN);
	int nz   = sm->Dim(SLC);

	// Full density k-spaces 
	Matrix<raw> tmp_in  (nx, ny, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); 
	Matrix<raw> tmp_sm  (nx, ny, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); 
	Matrix<raw> tmp_out (nx, ny, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);

	for (int i = 0; i < nz; i++ )
		for (int j = 0; j < nc; j++) {

			memcpy (&tmp_in.at(0), &in->at(0, 0, nc, 0, 0, 0, 0, 0, 0, nz, 0, 0, 0, 0, 0, 0),   nx*ny*sizeof(double)); 
			memcpy (&tmp_sm.at(0), &sm->at(0, 0, nc, 0, 0, 0, 0, 0, 0, nz, 0, 0, 0, 0, 0, 0),   nx*ny*sizeof(double));

			//((noncart::nufft*)ncs)->forward_2d(&tmp_in, &tmp_out);

			memcpy (&out->at(0, 0, nc, 0, 0, 0, 0, 0, 0, nz, 0, 0, 0, 0, 0, 0), &tmp_out.at(0), nx*ny*sizeof(double));

		}

	return OK;

}


/**
 * @brief               Compute right hand side (i.e. Multiply EH, Hermitian counterpart to E, with k-space data)
 *
 * @param  in           K-space samples along trajectory O (Nk x Nc)
 * @param  sm           Sensitivity maps                 O (Nx x Ny x Nc)
 * @param  kt           k-space trajectory               O (Nk x 1)
 * @param  ncs          Non-Cartesian strategy for non uniform ft
 * @param  out          Returned product                 O (Nx x Ny)
 */
RRSModule::error_code
EH (Matrix<raw>* in, Matrix<raw>* sm, Matrix<raw>* kt, solver_plan_complex* plan, double epsilon, int maxit, Matrix<raw>* out) {

	int ncoils   = sm->Dim(CHA);
	int nsamples = in->Size(); 

	double t=nfft_second();
	
	/* inverse trafo */
	solver_before_loop_complex(plan);
	
	for(int l=0; l < maxit; l++)  {
		/* break if dot_r_iter is smaller than epsilon*/
		if(plan->dot_r_iter<epsilon)
			break;
#ifdef VERBOSE
		fprintf(stderr,"%e,  %i of %i\n",sqrt(plan->dot_r_iter), l+1, maxit);
#endif
		solver_loop_one_step_complex(plan);
	}

	t=nfft_second()-t;

	/*
	  for (k=0;k<my_plan.N_total;k++) {
	  fprintf(fout_real,"%le ", creal(my_iplan.f_hat_iter[k]));
	  fprintf(fout_imag,"%le ", cimag(my_iplan.f_hat_iter[k]));
	  }
	*/
	
	return OK;
	
}

// L2 norm
bool Converged (Matrix<raw>* r, Matrix<raw>* a, double cgconv) {

	float rn;
	float an;

	r->norm(&rn);
	a->norm(&an);

	return((rn/an) < cgconv);

}

RRSModule::error_code
CGSENSE::Process () {
	
	Matrix<raw> a, b, p, q, r, r_new, k;
	
	double residue = 1;
	double delta   = 0;
	
	// CG iterations
	for (int iter = 0; iter < m_iter; iter++) {
		
		float       rn;
		float       an;
		float       rnewn;
		raw         rtmp;
		Matrix<raw> mtmp;

		r.norm(&rn);
		a.norm(&an);

		delta = rn/an;

		if (delta < m_cgconv)
			break;
		
		/* q     = eh(e(p , sensitivity, k), sensitivity, k); */
		E  (&p,    &m_sens, &k, &m_fplan, &mtmp);
		EH (&mtmp, &m_sens, &k, &m_iplan, m_epsilon, m_maxit, &q);

		/* b     = b + r(:)'*r(:)/(p(:)'*q(:))*p;*/
		rtmp  = (rn / (p.dotc(q)));
		mtmp  = p * rtmp;
		b     = b + mtmp;
		
		/* r_new = r - r(:)'*r(:)/(p(:)'*q(:))*q; */
		mtmp  = q * rtmp;
		r_new = r - mtmp;

		/* p     = r_new + r_new(:)'*r_new(:)/(r(:)'*r(:))*p */
		r_new.norm(&rnewn);
		rtmp  = rnewn/rn;
		mtmp  = p * rtmp;
		p     = r_new + mtmp;

		/*  r     = r_new; */ 
		r     = r_new;

	}
	
}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

