#include "SENSE.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

using namespace RRStrategy;

codeare::error_code 
SENSE::Init () {

	printf ("Intialising Cartesian SENSE ...\n");

	Attribute ("compgfm", &m_compgfm);
	Attribute ("nthreads",  &m_nthreads);
	Attribute ("lambda", &m_lambda);

	printf ("... done.\n\n");

	return codeare::OK;

}


codeare::error_code 
SENSE::Prepare () {

	printf ("Preparing Cartesian SENSE ...\n");

    AddMatrix<cxfl> ("image");

	printf ("  allocating Cartesian SENSE operator: ... "); fflush(stdout);

	Params p;

    p.Set("smaps", Get<cxfl>("sensitivities"));       // Sensitivities
    p.Set("fdims", size(Get<cxfl>("aliased"))); // Folded dimensions

	p.Set("compgfm", m_compgfm);              // Compute g-factors?
	p.Set("nthreads", m_nthreads);            // Number of threads for FFT
	p.Set("lambda", m_lambda);                // Tikhonov regularizing term

	m_cs = CSENSE<float> (p);                 // Cartesian SENSE operator

	printf ("... done.\n\n");

	return codeare::OK;

}


codeare::error_code
SENSE::Process () { 

    SimpleTimer st ("SENSE");

    Matrix<cxfl>& out = Get<cxfl>("image");
    Matrix<cxfl>& in = Get<cxfl>("aliased");

    out = m_cs ->* in;
    
    return codeare::OK;

}


codeare::error_code
SENSE::Finalise () { 
	return codeare::OK;
}



// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new SENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}



