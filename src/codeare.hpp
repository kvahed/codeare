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
#include <sstream>

char  *name = 0, *base_dir = 0, *config = 0, *debug = 0, *EMPTY = (char*)"", *ZERO = (char*)"0",
    *CURRENT = (char*)".", *data = 0, *infile = 0, *outfile = 0;
std::string config_file_uri;
bool   pulses;
int    rank;

bool commandline_opts (int argc, char** argv) {
  
#ifdef HAVE_MPI
  //Grid& gd = Grid::Instance();
  //if (gd.rk == 0) {
#endif
    cout << endl << "codeare " << VERSION ; 

#ifdef GIT_SHA1
    cout << " [" << GIT_SHA1 << "]";
#endif
    cout << endl;
#ifdef HAVE_MPI
//  }
#endif

    std::stringstream intrinsics;
#ifdef __MMX__
    intrinsics << "MMX ";
#endif
#ifdef __SSE__
    intrinsics << "SSE ";
#endif
#ifdef __SSE2__
    intrinsics << "SSE2 ";
#endif
#ifdef __SSSE3__
    intrinsics << "SSSE3 ";
#endif
#ifdef __AVX__
    intrinsics << "AVX ";
#endif
#ifdef __AVX2__
    intrinsics << "AVX2 ";
#endif
#ifdef __FMA__
    intrinsics << "FMA ";
#endif
#ifdef __FMA4__
    intrinsics << "FMA4 ";
#endif
    
    Options opt;
	
    opt.addUsage  ("Copyright (C) 2010-2014");
    opt.addUsage  ("Kaveh Vahedipour<kaveh@codeare.org>");
    opt.addUsage  ("NYU School of Medicine");
    opt.addUsage  (intrinsics.str().c_str());
    opt.addUsage  ("Usage:");
    opt.addUsage  ("codeare -c <configuration_file> -b <base_dir> [OPTIONS]");
    opt.addUsage  ("");
    opt.addUsage  (" -n, --name     Remote service name (optional)");
    opt.addUsage  (" -d, --debug    Debug level 0-40 (default: 0)");
    opt.addUsage  (" -b, --base     Base directory of files");
    opt.addUsage  ("                (default: current directory, i.e. '.'");
    opt.addUsage  (" -c, --config   Configuration XML file name (default: codeare.xml");
    opt.addUsage  (" -i, --infile   Binary data input file name");
    opt.addUsage  ("                (optional: May be assigned in configuration)");
    opt.addUsage  (" -o, --outfile  Binary data output file name");
    opt.addUsage  ("                (optional: May be assigned in configuration)");
    opt.addUsage  ("");
    opt.addUsage  (" -h, --help     Print this help screen");
    opt.addUsage  ("");
    
    opt.setFlag   ("help"  , 'h');
    opt.setOption ("name"  , 'n');
    opt.setOption ("debug" , 'd');
    opt.setOption ("config", 'c');
    opt.setOption ("base"  , 'b');
    
    opt.processCommandArgs(argc, argv);
    
    // Help screen
    if (/* !(opt.hasOptions()) ||*/ opt.getFlag("help")) {
        opt.printUsage();
        return false;
    } 
	
    // Debug level
    char* tmp = opt.getValue("debug");
    debug  = (tmp && atoi(tmp) >= 0 && atoi(tmp)  <= 40) ? tmp : ZERO;
    
    // Remote service's CORBA name default "" if specified, remote access is assumed
    tmp      = opt.getValue("name");
    name     = (tmp && strcmp(tmp,EMPTY)) ? tmp : 0;
    
    // Base directory for data
    tmp      = opt.getValue("base");
    base_dir = (tmp && strcmp(tmp,EMPTY)) ? tmp : CURRENT;
    
    // Configuration file, default: config.xml
    tmp      = opt.getValue("config");
    config   = (tmp && strcmp(tmp,EMPTY)) ? tmp : (char*) "config.xml";
    
    config_file_uri  = base_dir;
    config_file_uri += "/";
    config_file_uri += std::string(config);
    
    return true;
	
}

