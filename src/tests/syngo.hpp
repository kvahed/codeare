#include "matrix/IO.hpp"
#include "io/RawParser.hpp"

template <class T> bool
syngotest (Connector<T>* rc) {

	using namespace matrix::io::SyngoMR;
	
	std::string in = std::string (base + std::string (data));
	
	Matrix<cxfl>  meas;
	Matrix<float> sync;
	std::string   fname;
	
	if (rc->Attribute("filename", &fname) != TIXML_SUCCESS) {
		printf ("\n  No filename attribute specified. Exiting!\n");
		return false;
	}

	RawParser rp = RawParser (fname, true);
	if (rp.Status()!=MF_OK) {
		printf ("\n  Couldn't initialise RawParser. Exiting!\n");
		return false;
	}

	rp.Parse();
	
	Data<cxfl> dmeas = rp.GetMeas();
	meas = Matrix<cxfl> (dmeas.dims, dmeas.ress);
	meas.Dat() = dmeas.data;
	meas = squeeze(meas);

	printf ("    Matrix dims: %s\n", DimsToCString(meas));

	rp.CleanUp();

	return true;
	

}
