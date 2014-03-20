#include "IO.hpp"
#include "io/IOContext.hpp"

template <class T> bool
ismrmrdtest (Connector<T>* rc) {

	using namespace codeare::matrix::io;
	
	Params params;
	params.Set("scheme", std::string("/usr/local/ismrmrd/schema/ismrmrd.xsd"));
	IOContext ioc = IOContext (std::string (base_dir + std::string (data)), ISMRM, READ, params);

	//ioc.Read("dataset");

	return true;
	

}
