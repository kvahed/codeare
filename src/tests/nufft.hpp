#include "matrix/io/IOContext.hpp"

template <class T> bool
nuffttest (Connector<T>* rc) {

	using namespace codeare::matrix::io;

	Matrix<cxdb>   rawdata;
	Matrix<double> weights;
	Matrix<double> kspace;
	Matrix<cxdb>   img;

	std::string ipf = base;
	ipf += rc->GetElement("/config/data-in")->Attribute("fname");

	rc->Init(test);

	IOContext in = fopen (ipf, "r");
	rawdata = fread<cxfl> (in, "data");
	kspace  = fread<float> (in, "kspace");
	weights = fread<float> (in, "weights");
	fclose (in);
	
	rc->SetMatrix    ("data",    rawdata);
	rc->SetMatrix    ("weights", weights);
	rc->SetMatrix    ("kspace",  kspace);
	
	rc->Prepare      (test);
	rc->Process      (test);
	
	rc->GetMatrix    ("img", img);
	img.SetClassName("img");

	rc->Finalise(test);

	std::string opf = base;
	opf += rc->GetElement("/config/data-out")->Attribute("fname");
	IOContext out = fopen (opf, "w");
	fwrite (out, img);
	fclose (out);

	return true;
	
}
