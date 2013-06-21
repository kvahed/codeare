/**************
 ** includes **
 **************/
# include "Matrix.hpp"
# include "IOContext.hpp"
# include "DWT.hpp"

/**********************
 ** type definitions **
 **********************/
typedef double elem_type;


using namespace codeare::matrix::io;



/*******************
 ** test function **
 *******************/
int
main            (int argc, const char ** args)
{

    Configurable conf;
    conf.ReadConfig (args [1]);

	// Intro
	std::cout << std::endl;
	std::cout << " * Running test: \"dwt2\"!";
	std::cout << std::endl;
	std::cout << std::endl;


	// create oclMatrix from input file
	Matrix <elem_type> mat_in;
	IOContext ioc (rc->GetElement ("/config/data/in"), base, READ);
	mat_in = ioc.Read <elem_type> (rc->GetElement ("/config/data/in/signal"));

	// read DWT params
	wlfamily wl_fam = (wlfamily) atoi (rc->GetElement ("/config")->Attribute ("wl_fam"));
	const int wl_mem = atoi (rc->GetElement ("/config")->Attribute ("wl_mem"));

	int iterations = 5;

	// do something
	DWT <elem_type> dwt (mat_in.Dim (0), wl_fam, wl_mem, 4);
	Matrix <elem_type> mat_out_dwt = dwt * mat_in;
	Matrix <elem_type> mat_out_dwt_recon;

	for (int i = 0; i < iterations; i++)
	{
	    mat_out_dwt_recon = dwt ->* mat_out_dwt;
	    mat_out_dwt = dwt * mat_out_dwt_recon;
	}


	// output oclMatrix to output file

	IOContext ioc2 (rc->GetElement ("/config/data/out"), base, WRITE);
	ioc2.Write(mat_out_dwt, rc->GetElement ("/config/data/out/res-dwt"));
    ioc2.Write(mat_out_dwt_recon, rc->GetElement ("/config/data/out/res-dwt-recon"));

	// Outro
	std::cout << std::endl;
	std::cout << " * Finished test: \"dwt2\"!";
	std::cout << std::endl;
	std::cout << std::endl;

}
