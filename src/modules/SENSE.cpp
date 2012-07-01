#include "SENSE.hpp"

using namespace RRStrategy;

RRSModule::error_code 
SENSE::Init () {

	printf ("Intialising Cartesian SENSE ...\n");

	Attribute ("compgfm", &m_compgfm);
    printf ("  compute g-factor maps: %s \n", (m_compgfm) ? "true" : "false");
	Attribute ("ncpus",   &m_ncpus);
	printf ("  # threads:             %i \n", m_ncpus);

	printf ("... done.\n\n");

	return RRSModule::OK;

}


RRSModule::error_code 
SENSE::Prepare () {

	printf ("Preparing Cartesian SENSE ...\n");

	Matrix<cxfl>& smaps = GetCXFL("smaps");
	Matrix<cxfl>& fimgs = GetCXFL("fimgs");

	// We expect the first phase encoding dimension 
	// to be accelerated
	m_af = size(smaps, 1) / size(fimgs, 1);
	printf ("  # channels:            %zu \n", size(smaps,0));
	printf ("  # acceleration factor: %i \n", m_af);

	Matrix<cxfl>& image = AddMatrix 
		("image", (Ptr<Matrix<cxfl> >) NEW (Matrix<cxfl>(1)));

	printf ("  allocating Cartesian SENSE operator: ... "); fflush(stdout);
	m_cs = new CSENSE<float> (smaps, m_af, m_compgfm);
	printf ("done\n");
	printf ("... done.\n\n");

	return RRSModule::OK;

}


RRSModule::error_code
SENSE::Process () { 

	ticks cgstart = getticks();
	
	printf ("Processing SENSE ...\n");

	GetCXFL("image") = *m_cs ->* GetCXFL("fimgs");

	printf ("... done. WTime: %.4f seconds.\n\n", elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());

	return RRSModule::OK;

}


RRSModule::error_code
SENSE::Finalise () { 

	delete m_cs;

	return RRSModule::OK;

}



// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new SENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}



