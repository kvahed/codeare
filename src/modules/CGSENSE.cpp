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

void 
CGSENSE::E  (Matrix< raw >* in, Matrix< raw >* out) {
	
	int ncoils   = m_sens.Dim(CHA);
	int nsamples = m_raw.Size(); 

	// Full density k-spaces per coil
	Matrix<raw> fdkspaces;
	out->Reset(m_raw.Dim());
	
	// Reverse grid full density k-space on actual trajectory
	

}

void 
CGSENSE::EH (Matrix< raw >* in, Matrix< raw >* out) {

	int ncoils   = m_sens.Dim(CHA);
	int nsamples = m_raw.Size(); 

	// Backward nufft

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

