#include "IO.hpp"
#include "io/IOContext.hpp"

template <class T> bool
syngotest (Connector<T>* rc) {

	/*
	struct Params {

		short int mode;

		bool      undo_os;
		bool      rep_wise;
		bool      average;
		bool      ext_ref_only;

		Params () {
			mode         = TSE;
			undo_os      = 1;
			rep_wise     = 0;
			average      = 1;
			ext_ref_only = 0;
		}

	};
*/



	using namespace codeare::matrix::io;

	Params params;
	params.Set("undo_os", true);
	params.Set("res_wise", false);
	params.Set("average", (int)1);
	params.Set("ext_ref_only", false);
	params.Set("mode", std::string("TSE"));

	std::string in = std::string (base + std::string (data));
	
	Matrix<cxfl> meas;
	Matrix<float> sync;
	std::string   fname;
	
	if (rc->Attribute("filename", &fname) != TIXML_SUCCESS) {
		printf ("\n  No filename attribute specified. Exiting!\n");
		return false;
	}

	IOContext ioc (fname, SYNGOMR, params);

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
