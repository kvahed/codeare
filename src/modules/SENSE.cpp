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

	AddMatrix ("image", (Ptr<Matrix<cxfl> >) NEW (Matrix<cxfl>()));

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

	ticks cgstart = getticks();
	
	printf ("Processing SENSE ...\n");

	Get<cxfl>("image") = *m_cs ->* Get<cxfl>("fimgs");

	printf ("... done. WTime: %.4f seconds.\n\n",
			elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());

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



