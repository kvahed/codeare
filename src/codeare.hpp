#include "options.h"
#include "GitSHA1.hpp"

#include "Connector.hpp"

#ifndef __WIN32__
#    include "config.h"
#endif

using namespace std;
using namespace RRClient;

#ifdef HAVE_MPI
#  include "Grid.hpp"
#endif

#include <time.h>
#include <stdio.h>

char  *name = 0, *base_dir = 0, *config = 0, *verbose = 0, *EMPTY = (char*)"", *ZERO = (char*)"0";
const char *test;
std::string config_file_uri;
bool   pulses;
int    rank;

bool commandline_opts (int argc, char** argv) {

#ifdef HAVE_MPI
  //Grid& gd = Grid::Instance();
  //if (gd.rk == 0) {
#endif
	cout << "codeare " << VERSION  << endl;

#ifdef GIT_COMMIT
	cout << "Commit " << g_GIT_SHA1 << " [" << GIT_COMMIT_DATE << "]" << endl;
#endif
#ifdef HAVE_MPI
//  }
#endif

	Options opt;
		
	opt.addUsage  ("Copyright (C) 2010-2014");
	opt.addUsage  ("Kaveh Vahedipour<k.vahedipour@fz-juelich.de>");
	opt.addUsage  ("Juelich Research Centre");
	opt.addUsage  ("Medical Imaging Physics");
	opt.addUsage  ("");
	opt.addUsage  ("Usage:");
	opt.addUsage  ("codeare -c <configuration_file> -b <base_dir> [OPTIONS]");
	opt.addUsage  ("");
	opt.addUsage  (" -n, --name    Remote service name (for example: ReconService)");
	opt.addUsage  (" -v, --verbose Debug level 0-40 (default: 0)");
	opt.addUsage  (" -b, --base    Base directory of approved files.");
	opt.addUsage  (" -c, --config  Configuration XML (NuFFT, CGSENSE).");
	opt.addUsage  ("");
	opt.addUsage  (" -h, --help    Print this help screen");
	opt.addUsage  ("");
	
	opt.setFlag   ("help"    ,'h');
	opt.setOption ("name"   , 'n');
	opt.setOption ("verbose", 'v');
	opt.setOption ("config" , 'c');
	opt.setOption ("base"   , 'b');
	
	opt.processCommandArgs(argc, argv);
	
	// Help screen
	if ( !(opt.hasOptions()) || opt.getFlag("help")) {
		opt.printUsage();
		return false;
	} 
	
	// Debug level
	verbose = (opt.getValue("verbose" )         &&
		  atoi(opt.getValue("verbose"))  >= 0   &&
		  atoi(opt.getValue("verbose"))  <= 40) ? opt.getValue("verbose") : ZERO;

	// Remote service's CORBA name default ""
	name    = (opt.getValue("name")   &&
			   opt.getValue("name")   != EMPTY) ? opt.getValue(   "name") : 0;

	// Base directory for data
	base_dir    = (opt.getValue("base")   &&
			   opt.getValue("base")   != EMPTY) ? opt.getValue(   "base") : 0;

	// Configuration file
	config  = (opt.getValue("config") &&
			   opt.getValue("config") != EMPTY) ? opt.getValue( "config") : 0;

	config_file_uri  = base_dir;
	config_file_uri += "/";
	config_file_uri += std::string(config);

	return true;
	
}


/*#ifdef HAVE_MAT_H
bool init (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//std::string("codeare").c_str();
	return true;
}
#endif*/
