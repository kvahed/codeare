/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
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

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <signal.h>
                                                                                
#include "ReconServant.h"

#ifndef __WIN32__
#include "config.h"
#endif

#include "options.h"

#ifndef SVN_REVISION
	#define SVN_REVISION "unkown"
#endif

char*  name;
char*  debug;
char*  logfile;

using namespace std;
using namespace RRServer;

bool init (int argc, char** argv);
                                                                                
int main (int argc, char** argv) {

	if (!init (argc, argv))
		return 1;

	// --------------------------------------------------------------------------
	// Start CORBA server:
	// --------------------------------------------------------------------------
	try {

		streambuf* out;
		streambuf* err;
		ofstream   log (logfile);

		if (log.is_open()) {
			 out = cout.rdbuf(log.rdbuf());
			 err = cerr.rdbuf(log.rdbuf());
		} else {
			cout << "Could not open logfile " << logfile << "." << endl;
			cout << "Exiting :(" << endl << endl;
			return 1;
		}
		
		// Initialise ORB
		const char*    options[][2] = { { (char*)"traceLevel", debug }, { 0, 0 } };
		CORBA::ORB_var orb          = CORBA::ORB_init(argc, argv, "omniORB4", options);
		
		// Get reference to the RootPOA.
		CORBA::Object_var obj = orb->resolve_initial_references("RootPOA");
		PortableServer::POA_var _poa = PortableServer::POA::_narrow(obj.in());

		// Initialise servant
		ReconServant* myReconServant = new ReconServant();

		// Activate in RootPDA
		PortableServer::ObjectId_var myReconServant_oid
			= _poa->activate_object(myReconServant);
                                                                                
		// Obtain object reference from servant and register in naming service
		CORBA::Object_var SA_obj = myReconServant->_this();
                                                                                
		// Obtain a reference to the object, and print it out as string IOR.
		CORBA::String_var sior(orb->object_to_string(SA_obj.in()));
                                                                                
		// Bind to the name server and lookup 
		CORBA::Object_var obj1=orb->resolve_initial_references("NameService");
		assert(!CORBA::is_nil(obj1.in()));

		// Get context
		CosNaming::NamingContext_var nc = CosNaming::NamingContext::_narrow(obj1.in());
		assert(!CORBA::is_nil(nc.in()));
                                                                                
		// Resolve name
		CosNaming::Name m_name;
		m_name.length(1);
		m_name[0].id=CORBA::string_dup(name);
		nc->rebind (m_name,SA_obj.in());
                                                                                
        // Activate POA manager
		PortableServer::POAManager_var pmgr = _poa->the_POAManager();
		pmgr->activate();
                                                                                
		// Accept requests from clients
		orb->run();
                                                                                
		// If orb leaves event handling loop.
		// - currently configured never to time out (??)
		orb->destroy();
                                                                                
		free(m_name[0].id); // str_dup does a malloc internally

		cout.rdbuf(out);
		cerr.rdbuf(err);

	}
	
	catch(CORBA::SystemException&) {
		cerr << "Caught CORBA::SystemException." << endl;
	}
	catch(CORBA::Exception&) {
		cerr << "Caught CORBA::Exception." << endl;
	}
	catch(omniORB::fatalException& fe) {
		cerr << "Caught omniORB::fatalException:" << endl;
		cerr << "  file: " << fe.file() << endl;
		cerr << "  line: " << fe.line() << endl;
		cerr << "  mesg: " << fe.errmsg() << endl;
	}
	catch(...) {
		cerr << "Caught unknown exception." << endl;
	}
	
	return 0;

}


bool init (int argc, char** argv) {

	cout << endl;
#ifdef VERSION
	cout << "jrrs "         << VERSION                                        << endl;
#else
	cout << "jrrs "         << endl;
#endif
	cout << "juelich remote reconstruction service "                          << endl;
	cout << "Reconstruction server "  << " [build " << SVN_REVISION << "]"    << endl;
    cout << "Copyright (C) 2010	Kaveh Vahedipour<k.vahedipour@fz-jeulich.de>" << endl;
	cout << "Kaveh Vahedipour -  k.vahedipour@fz-juelich.de"                  << endl;
	cout << "Juelich Research Centre"                                         << endl;
	cout << "Institute of Neuroscience and Medicine"                          << endl;
	cout << "Medical Imaging Physics"                                         << endl;
	cout << endl;

	Options *opt = new Options();
	
	opt->addUsage  ("Usage: testclt --name <name> [OPTIONS]");
	opt->addUsage  ("");
	opt->addUsage  (" -n, --name    Remote service name (for example RemoteRecon)");
	opt->addUsage  (" -d, --debug   Debug level 0-40 (default: 5)");
	opt->addUsage  (" -l, --logfile Log file (default: ./reconserver.log)");
	opt->addUsage  ("");
	opt->addUsage  (" -h, --help    Print this help screen"); 
	opt->addUsage  ("");
	
	opt->setFlag   ("help"    , 'h');

	opt->setOption ("logfile" , 'l');
	opt->setOption ("debug"   , 'd');
	opt->setOption ("name"    , 'n');
	
	opt->processCommandArgs(argc, argv);
	
	if ( !(opt->hasOptions()) || opt->getFlag("help") ) {
		
		opt->printUsage();
		delete opt;
		return true;
		
	} 
	
	debug   = (opt->getValue("debug")                && 
			   atoi(opt->getValue("debug")) >= 0     &&  
			   atoi(opt->getValue("debug")) <= 40)    ? opt->getValue("debug")   : (char*)"5";

	logfile = (opt->getValue("logfile")              &&
			   opt->getValue("logfile") != (char*)"") ? opt->getValue("logfile") : (char*)"./reconserver.log";

	name    = (opt->getValue("name")                 &&      
			   opt->getValue("name") != (char*)"")    ? opt->getValue("name")    : (char*)"ReconService";

	delete opt;
	return true;

}
