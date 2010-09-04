#include "CGSENSE.h"

using namespace RRStrategy;

CGSENSE::CGSENSE () {

	// Copy incoming data to m_temp             */
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
	
	N[0]  = m_raw.Dim(COL);
	N[1]  = m_raw.Dim(LIN);
	N[2]  = m_raw.Dim(SLC);

	Nk[0] = m_helper.Dim(COL);
	Nk[1] = m_helper.Dim(LIN);
	Nk[2] = m_helper.Dim(SLC);

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
	
	//m_nufft.init (2, N, Nk, m_helper.at(0), 2, 2.4, 32, 0, 0);


}


/**
 * @brief               Compute left hand side (i.e. Multiply E with spatial (image) data)
 *
 * @param  in           Original discretised sample O (Nx x Ny)
 * @param  sm           Sensitivity maps            O (Nx x Ny x Nc)
 * @param  kt           k-space trajectory          O (Nk x 1)
 * @param  out          Result                      O (Nk x Nc)
 */
RRSModule::error_code 
E  (Matrix<raw>* in, Matrix<raw>* sm, Matrix<raw>* kt, noncart::strategy* ncs, Matrix<raw>* out) {
	
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

			((noncart::nufft*)ncs)->forward_2d(&tmp_in, &tmp_out);

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
 * @param  out          Returned product                 O (Nx x Ny)
 */
RRSModule::error_code
EH (Matrix<raw>* in, Matrix<raw>* sm, Matrix<raw>* kt, noncart::strategy* ncs, Matrix<raw>* out) {

	int ncoils   = sm->Dim(CHA);
	int nsamples = in->Size(); 

	return OK;

}

RRSModule::error_code
CGSENSE::Process () {

	Matrix<raw> *p, *s, *k, *q;

	// CG iterations
	for (int iter = 0; iter < m_iter; iter++) {

		EH (p, s, k, &m_nufft, q);

		delete p;
		p = new Matrix<raw> (*(q));

		E  (p, s, k, &m_nufft, q);

		/*
		  delta = r(:)'*r(:)/(a(:)'*a(:));
		  q     = eh(e(p , sensitivity, k), sensitivity, k);
		  b     = b + r(:)'*r(:)/(p(:)'*q(:))*p;
		  r_new = r - r(:)'*r(:)/(p(:)'*q(:))*q;
		  p     = r_new + r_new(:)'*r_new(:)/(r(:)'*r(:))*p;
		  r     = r_new;
		*/
		
	}

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

