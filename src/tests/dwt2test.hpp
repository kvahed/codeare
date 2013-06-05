/**************
 ** includes **
 **************/
# include "matrix/Matrix.hpp"
# include "matrix/io/IOContext.hpp"
# include "matrix/dwt/DWT.hpp"
# include "matrix/dwt/DWT2.hpp"

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
	std::cout << " * Running test: \"dwt2\"!";
	std::cout << std::endl;
	std::cout << std::endl;


	// create oclMatrix from input file
	Matrix <elem_type> mat_in;
	IOContext ioc (rc->GetElement ("/config/data/in"), base, READ);
	mat_in = ioc.Read <elem_type> (rc->GetElement ("/config/data/in/signal"));


	// do something
	DWT <elem_type> dwt (mat_in.Dim (0), WL_DAUBECHIES, 4);
	Matrix <elem_type> mat_out_dwt = dwt * mat_in;

	DWT2 <elem_type> dwt2 (2);
	Matrix <elem_type> mat_out_dwt2 = dwt2 * mat_in;
	Matrix <elem_type> mat_out_dwt2_recon = dwt2 ->* mat_out_dwt2;


	// output oclMatrix to output file

//	std::string ofname (base + std::string (rc->GetElement ("/config/data/out") -> Attribute ("fname")));
//	std::string odname (rc->GetElement ("/config/data/out") -> Attribute ("dname"));

	IOContext ioc2 (rc->GetElement ("/config/data/out"), base, WRITE);
	ioc2.Write(mat_out_dwt, rc->GetElement ("/config/data/out/res-dwt"));
	ioc2.Write(mat_out_dwt2, rc->GetElement ("/config/data/out/res-dwt2"));
	ioc2.Write(mat_out_dwt2_recon, rc->GetElement ("/config/data/out/res-dwt2-recon"));

	// Outro
	std::cout << std::endl;
	std::cout << " * Finished test: \"dwt2\"!";
	std::cout << std::endl;
	std::cout << std::endl;

}
