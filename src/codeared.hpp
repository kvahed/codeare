/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301  USA
 */


#ifndef CODEARED_H_
#define CODEARED_H_


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <signal.h>



#include "ReconServant.hpp"
#include "CorbaExceptions.hpp"

#ifndef __WIN32__
#include "config.h"
#endif

#include "options.h"
#include "GitSHA1.hpp"

#ifndef SVN_REVISION
    #define SVN_REVISION "unkown"
#endif

char  *name, *debug, *logfile, *port, *EMPTY = (char*)"", *FIVE = (char*)"5";

using namespace std;
using namespace RRServer;

bool init (int argc, char** argv) {

#ifdef HAVE_MPI
  //Grid& gd = Grid::Instance();
  //if (gd.rk == 0) {
#endif
    cout << endl << "codeare service " << VERSION ; 

#ifdef GIT_SHA1
    cout << " [" << GIT_SHA1 << "]";
#endif
    cout << endl;
#ifdef HAVE_MPI
//  }
#endif

    Options *opt = new Options();

    opt->addUsage  ("Copyright (C) 2010-2014");
    opt->addUsage  ("Kaveh Vahedipour<k.vahedipour@fz-juelich.de>");
    opt->addUsage  ("Juelich Research Centre");
    opt->addUsage  ("Medical Imaging Physics");
    opt->addUsage  ("");
    opt->addUsage  ("Usage:");
    opt->addUsage  ("codeared -n <servicename> [-l <logfile> -d <debuglevel>]");
    opt->addUsage  ("");
    opt->addUsage  (" -n, --name     Service name (default: codeared)");
    opt->addUsage  (" -d, --debug    Debug level 0-40 (default: 5)");
    opt->addUsage  (" -l, --logfile  Log file (default: ./reconserver.log)");
    opt->addUsage  (" -p, --httpport http service port (default 8080)");
    opt->addUsage  ("");
    opt->addUsage  (" -h, --help     Print this help screen");

    opt->setFlag   ("help"    , 'h');

    opt->setOption ("logfile" , 'l');
    opt->setOption ("debug"   , 'd');
    opt->setOption ("name"    , 'n');
    opt->setOption ("httpport", 'p');


    opt->processCommandArgs(argc, argv);

    if (/* !(opt->hasOptions()) ||*/ opt->getFlag("help") ) {

        opt->printUsage();
        delete opt;
        return false;

    }

    char* tmp = opt->getValue("debug");
    debug   = (tmp && atoi(tmp) >= 0 && atoi(tmp)  <= 40) ? tmp : FIVE;

    tmp = opt->getValue("logfile");
    logfile = (tmp && strcmp(tmp,EMPTY)) ? tmp : (char*)"./codeared.log";

    tmp = opt->getValue("name");
    name    = (tmp && strcmp(tmp,EMPTY)) ? tmp : (char*) "codeared";

    tmp = opt->getValue("httpport");
    port    = (tmp && atoi(tmp) >= 0 && atoi(tmp) <= 65536) ? tmp : (char*)"8080";

    delete opt;
    return true;

}




#endif /* CODEARED_H_ */
