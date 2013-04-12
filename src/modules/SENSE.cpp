#include "SENSE.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

using namespace RRStrategy;

error_code 
SENSE::Init () {

	printf ("Intialising Cartesian SENSE ...\n");

	Attribute ("compgfm", &m_compgfm);
    printf ("  compute g-factor maps: %s \n", (m_compgfm) ? "true" : "false");
	Attribute ("ncpus",   &m_ncpus);
	printf ("  # threads:             %i \n", m_ncpus);

	printf ("... done.\n\n");

	return OK;

}


error_code 
SENSE::Prepare () {

	printf ("Preparing Cartesian SENSE ...\n");

	Matrix<cxfl>& image = AddMatrix 
		("image", (Ptr<Matrix<cxfl> >) NEW (Matrix<cxfl>(1)));

	printf ("  allocating Cartesian SENSE operator: ... "); fflush(stdout);
	//m_cs = new CSENSE<float> (smaps, m_af, m_compgfm);
	printf ("done\n");
	printf ("... done.\n\n");

	Params p;
	p.Set("smaps_name", std::string("smaps"));
	p.Set("fimgs_name", std::string("fimgs"));
	p.Set("compgfm", m_compgfm);
	m_cs = new CSENSE<float> (p);

	return OK;

}


error_code
SENSE::Process () { 

	ticks cgstart = getticks();
	
	printf ("Processing SENSE ...\n");

	Get<cxfl>("image") = *m_cs ->* Get<cxfl>("fimgs");

	printf ("... done. WTime: %.4f seconds.\n\n", elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());

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



