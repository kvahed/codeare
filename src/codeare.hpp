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
	opt.addUsage  (" -b, --base    Base directory of approved files. Default: current directory ('.').");
	opt.addUsage  (" -c, --config  Configuration XML. Default: codeare.xml.");
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
	if (/* !(opt.hasOptions()) ||*/ opt.getFlag("help")) {
	    opt.printUsage();
	    return false;
	} 
	
	// Debug level
	char* tmp = opt.getValue("verbose");
	verbose  = (tmp && atoi(tmp) >= 0 && atoi(tmp)  <= 40) ? tmp : ZERO;
	
	// Remote service's CORBA name default "" if specified, remote access is assumed
	tmp      = opt.getValue("name");
	name     = (tmp && strcmp(tmp,EMPTY)==0) ? tmp : 0;
	
	// Base directory for data
	tmp      = opt.getValue("base");
	base_dir = (tmp && strcmp(tmp,EMPTY)==0) ? tmp : EMPTY;
	
	// Configuration file, default: config.xml
	tmp      = opt.getValue("config");
	config   = (tmp && strcmp(tmp,EMPTY)==0) ? tmp : (char*) "config.xml";
	
	config_file_uri  = base_dir;
	config_file_uri += "/";
	config_file_uri += std::string(config);

	return true;
	
}

