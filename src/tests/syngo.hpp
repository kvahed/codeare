#include "IO.hpp"
#include "io/IOContext.hpp"

template <class T> bool
syngotest (Connector<T>* rc) {

	using namespace codeare::matrix::io;
	
	std::string in = std::string (base + std::string (data));
	
	Matrix<cxfl> meas;
	Matrix<float> sync;
	std::string   fname;
	
	if (rc->Attribute("filename", &fname) != TIXML_SUCCESS) {
		printf ("\n  No filename attribute specified. Exiting!\n");
		return false;
	}

	IOContext ioc (fname, SYNGOMR);

	if (ioc.Status() != RRSModule::OK) {
		printf ("\n  Couldn't initialise RawParser. Exiting!\n");
		return false;
	}

	ioc.Read("");
	meas = squeeze(ioc.Get<cxfl>("meas"));
	MXDump (meas, "meas.mat", "meas");

	printf ("    Matrix dims: %s\n", DimsToCString(meas));

	return true;
	

}
