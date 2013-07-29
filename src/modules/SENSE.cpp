#include "SENSE.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

using namespace RRStrategy;

error_code 
SENSE::Init () {

	printf ("Intialising Cartesian SENSE ...\n");

	Attribute ("compgfm", &m_compgfm);
	Attribute ("nthreads",  &m_nthreads);
	Attribute ("lambda", &m_lambda);

	printf ("... done.\n\n");

	return OK;

}


error_code 
SENSE::Prepare () {

	printf ("Preparing Cartesian SENSE ...\n");

    AddMatrix<cxfl> ("image");

	printf ("  allocating Cartesian SENSE operator: ... "); fflush(stdout);

	Params p;

	p.Set("smaps_name", std::string("smaps"));
	p.Set("fimgs_name", std::string("fimgs"));
	p.Set("compgfm", m_compgfm);
	p.Set("nthreads", m_nthreads);
	p.Set("lambda",   m_lambda);

	m_cs = new CSENSE<float> (p);

	printf ("... done.\n\n");

	return OK;

}


error_code
SENSE::Process () { 

    SimpleTimer st ("SENSE");

    Matrix<cxfl>& out = Get<cxfl>("image");
    Matrix<cxfl>& in = Get<cxfl>("fimgs");
    CSENSE<float>& sense = *m_cs;

    out = sense ->* in;
    
	return OK;

}


error_code
SENSE::Finalise () { 

	delete m_cs;

	return OK;

}



// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new SENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}



