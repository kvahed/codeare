#include "CGSENSE.h"

CGSENSE::CGSENSE () {

	// Copy incoming data to m_temp
	m_temp (m_raw);

	// Reshape for outgoing format
	Attribute("Nx", &m_raw.Dim(0));
	Attribute("Ny", &m_raw.Dim(0));
	Attribute("Nz", &m_raw.Dim(0));

}

void 
CGSENSE::E  (Matrix* in, Matrix* out) {
	
	int ncoils   = m_sens.Dim(CHA);
	int nsamples = m_raw.Size(); 

	// Reset output to zeros
	out->Reset();
	

	return 

}

void 
CGSENSE::EH (Matrix* in, Matrix* out) {

	int ncoils   = m_sens.Dim(CHA);
	int nsamples = m_raw.Size(); 

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

