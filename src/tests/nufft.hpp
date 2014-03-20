#include "matrix/io/IOContext.hpp"

bool nuffttest (Connector& rc) {


	using namespace codeare::matrix::io;

	Matrix<cxdb>   rawdata;
	Matrix<double> weights;
	Matrix<double> kspace;
	Matrix<cxdb>   img;

	std::string ipf = base_dir;
	ipf += rc.GetElement("/config/data-in")->Attribute("fname");

	rc.Init(test);

	IOContext in = fopen (ipf);
	rawdata = fread<cxfl> (in, "data");
	kspace  = fread<float> (in, "kspace");
	weights = fread<float> (in, "weights");
	fclose (in);
	
	rc.SetMatrix    ("data",    rawdata);
	rc.SetMatrix    ("weights", weights);
	rc.SetMatrix    ("kspace",  kspace);
	
	rc.Prepare      (test);
	rc.Process      (test);
	
	rc.GetMatrix    ("img", img);

	rc.Finalise(test);

	std::string opf = base_dir;
	opf += rc.GetElement("/config/data-out")->Attribute("fname");

	IOContext out = fopen (opf, WRITE);
	fwrite (out, img);
	fclose (out);

	return true;
	
}
