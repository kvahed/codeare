#include "CGSENSE.h"

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

	

}


E  (Matrix<raw>* data_in, matrix<raw>* sensitivity, Matrix<raw>* k, Matrix<raw>* data_out) {
	
	int ncoils   = sensitivity->Dim(CHA);
	int nsamples = data_in->Size(); 

	// Full density k-spaces 
	Matrix<raw> FT;
	out->Reset(m_raw.Dim());

	// Reverse grid full density k-space on actual trajectory

	return 0;

}



EH (Matrix<raw>* data_in, matrix<raw>* sensitivity, Matrix<raw>* k, Matrix<raw>* data_out) {

	int ncoils   = sensitivity->Dim(CHA);
	int nsamples = data_in->Size(); 

}

RRSModules::error_code
Process () {

	// CG iterations
	for (int iter = 0; i < m_iter; i++) {

		EH (p, s, k, q);
		E  (q, s, k, q);
		
	}

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

