/**************
 ** includes **
 **************/
# include "Matrix.hpp"
# include "io/IOContext.hpp"
# include "DWT.hpp"
# include "Configurable.hpp"
# include <sstream>

/**********************
 ** type definitions **
 **********************/
typedef double elem_type;


using namespace codeare::matrix::io;



/*******************
 ** test function **
 *******************/
int
main            (int argc, char ** args)
{

    if (argc != 3)
    {
        std::cerr << " usage: t_dwt2 <base_dir> <config_file> " << std::endl;
        exit (-1);
    }

    char * base = args [1];

    std::stringstream ss;
    ss << base << args [2];

    Configurable conf;
    conf.ReadConfig (ss.str ().c_str ());

	// Intro
	std::cout << std::endl;
	std::cout << " * Running test: \"dwt2\"!";
	std::cout << std::endl;
	std::cout << std::endl;


	// create oclMatrix from input file
	Matrix <elem_type> mat_in;
	IOContext ioc (conf.GetElement ("/config/data/in"), base, READ);
	mat_in = ioc.Read <elem_type> (conf.GetElement ("/config/data/in/signal"));

	// read DWT params
	wlfamily wl_fam = (wlfamily) atoi (conf.GetElement ("/config")->Attribute ("wl_fam"));
	const int wl_mem = atoi (conf.GetElement ("/config")->Attribute ("wl_mem"));
	const int wl_scale = atoi (conf.GetElement ("/config")->Attribute ("wl_scale"));

	int iterations = 5;

	// do something
	DWT <elem_type> dwt (mat_in.Dim (0), wl_fam, wl_mem, wl_scale);
	Matrix <elem_type> mat_out_dwt = dwt * mat_in;
	Matrix <elem_type> mat_out_dwt_recon = dwt ->* mat_out_dwt;

	for (int i = 0; i < iterations; i++)
	{
	    mat_out_dwt = dwt * mat_out_dwt_recon;
        mat_out_dwt_recon = dwt ->* mat_out_dwt;
	}

	// output oclMatrix to output file

	IOContext ioc2 (conf.GetElement ("/config/data/out"), base, WRITE);
	ioc2.Write(mat_out_dwt, conf.GetElement ("/config/data/out/res-dwt"));
    ioc2.Write(mat_out_dwt_recon, conf.GetElement ("/config/data/out/res-dwt-recon"));

	// Outro
	std::cout << std::endl;
	std::cout << " * Finished test: \"dwt2\"!";
	std::cout << std::endl;
	std::cout << std::endl;

}
