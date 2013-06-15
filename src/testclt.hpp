#include "options.h"

#include "LocalConnector.hpp"
#if defined HAVE_OMNIORB4_CORBA_H && defined REMOTE
    #include "RemoteConnector.hpp"
#endif	

#ifndef __WIN32__
    #include "config.h"
#endif

#ifndef SVN_REVISION
	#define SVN_REVISION "unkown"
#endif

using namespace std;
using namespace RRClient;

#include <time.h>
#include <stdio.h>

char*  name;
char*  base;
char*  data;
char*  config;
char*  verbose;
char*  test;
bool   pulses;
int    rank;

#include "tests/tests.hpp"

bool init (int argc, char** argv) {

#ifdef HAVE_MPI
    Grid& gd = *Grid::Instance();
#endif

//	if (gd.rk == 0) {
		cout << endl;
#ifdef REMOTE
		cout << "codeare remote client " << VERSION  << endl;
#else
		cout << "codeare local client "  << VERSION  << endl;
#endif
        
#ifdef GIT_COMMIT
		cout << "Commit " << GIT_COMMIT << " [" << GIT_COMMIT_DATE << "]" << endl;
#endif
        //  }
    
		Options *opt = new Options();
		
		opt->addUsage  ("Copyright (C) 2010-2012");
		opt->addUsage  ("Kaveh Vahedipour<k.vahedipour@fz-juelich.de>");
		opt->addUsage  ("Juelich Research Centre");
		opt->addUsage  ("Medical Imaging Physics");
		opt->addUsage  ("");
		opt->addUsage  ("Usage:");
#ifdef REMOTE
		opt->addUsage  ("rclient -n <servicename> -t <test> [OPTIONS]");
#else
		opt->addUsage  ("lclient -t <test> [OPTIONS]");
#endif
		opt->addUsage  ("");
#ifdef REMOTE
		opt->addUsage  (" -n, --name    Remote service name (for example: ReconService)");
#endif
		opt->addUsage  (" -t, --test    Test case (default: DummyRecon. Just connectivity test)");
		opt->addUsage  (" -v, --verbose Debug level 0-40 (default: 0)");
		opt->addUsage  (" -b, --base    Base directory of approved files.");
		opt->addUsage  (" -c, --config  Configuration XML (NuFFT, CGSENSE).");
		opt->addUsage  (" -d, --data    Incoming binary data in HDF5 format.");
		opt->addUsage  (" -p, --pulses  Pulses (Only excitation).");
		opt->addUsage  ("");
		opt->addUsage  (" -h, --help    Print this help screen"); 
		opt->addUsage  ("");
		
		opt->setFlag   ("help"    ,'h');
		opt->setFlag   ("pulses"  ,'p');
#ifdef REMOTE
		opt->setOption ("name"   , 'n');
#endif
		opt->setOption ("test"   , 't');
		opt->setOption ("verbose", 'v');
		opt->setOption ("config" , 'c');
		opt->setOption ("data"   , 'd');
		opt->setOption ("base"   , 'b');
		
		opt->processCommandArgs(argc, argv);
		
		if ( !(opt->hasOptions()) || opt->getFlag("help")) {
			
			opt->printUsage();
			delete opt;
			return false;
			
		} 
		
		verbose = (opt->getValue("verbose" ) && 
				   atoi(opt->getValue("verbose")) >= 0 && 
				   atoi(opt->getValue("verbose")) <= 40) ? opt->getValue("verbose") : (char*) "0";
		name    = (opt->getValue("name")   && 
				   opt->getValue("name")   != (char*)"") ? opt->getValue(   "name") : (char*) "codeare" ;
		base    = (opt->getValue("base")   && 
				   opt->getValue("base")   != (char*)"") ? opt->getValue(   "base") : (char*) "";
		data    = (opt->getValue("data")   && 
				   opt->getValue("data")   != (char*)"") ? opt->getValue(   "data") : (char*) "";
		config  = (opt->getValue("config") && 
				   opt->getValue("config") != (char*)"") ? opt->getValue( "config") : (char*) "";
		test    = (opt->getValue("test"  ) && 
				   opt->getValue("test"  ) != (char*)"") ? opt->getValue(   "test") : (char*) "DummyRecon";
		
		delete opt;

		
	return true;
	
}


#ifdef HAVE_MAT_H

bool init (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	//std::string("codeare").c_str();
	return true;

}


#endif
