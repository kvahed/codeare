/**************
 ** includes **
 **************/
# include "matrix/Matrix.hpp"
# include "matrix/io/IOContext.hpp"
# include "matrix/dwt/DWT.hpp"

/**********************
 ** type definitions **
 **********************/
typedef double elem_type;


/*******************
 ** test function **
 *******************/
template <class T> bool
dwt2test (RRClient::Connector<T>* rc)
{

	// Intro
	std::cout << std::endl;
	std::cout << " * Running test: \"ocl_dwt\"!";
	std::cout << std::endl;
	std::cout << std::endl;


	// create oclMatrix from input file
	Matrix <elem_type> mat_in;
	IOContext ioc (rc->GetElement ("/config/data/in"), base, READ);
	std::cout << " bla " << std::endl;
	mat_in = ioc.Read <elem_type> (rc->GetElement ("/config/data/in"));


	// do something
	DWT <elem_type> dwt (WL_HAAR);
	mat_in = dwt * mat_in;


	// output oclMatrix to output file

//	std::string ofname (base + std::string (rc->GetElement ("/config/data/out") -> Attribute ("fname")));
//	std::string odname (rc->GetElement ("/config/data/out") -> Attribute ("dname"));

	IOContext ioc2 (rc->GetElement ("/config/data/out"), base, WRITE);
	ioc2.Write(mat_in, rc->GetElement ("/config/data/out"));

	// Outro
	std::cout << std::endl;
	std::cout << " * Finished test: \"ocl_dwt\"!";
	std::cout << std::endl;
	std::cout << std::endl;

}
