/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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

#include "testclt.hpp"

int main (int argc, char** argv) {
	
	if (init (argc, argv)) {

		Connector<RemoteConnector>* rc; 
		Connector<LocalConnector>*  lc;

		if (remote) 
			rc = new Connector<RemoteConnector> (name, verbose);
		else 
			lc = new Connector<LocalConnector>  (name, verbose);
		
		if      (!strcmp (test, "CGSENSE")              ) (remote) ?  cgsensetest (rc) :  cgsensetest (lc);
		else if (!strcmp (test, "DirectMethod")         ) (remote) ?       dmtest (rc) :       dmtest (lc);
		else if (!strcmp (test, "NuFFT")                ) (remote) ?    nuffttest (rc) :    nuffttest (lc); 
		else if (!strcmp (test, "NuFFT_OMP")            ) (remote) ?    nuffttest (rc) :    nuffttest (lc);
		else if (!strcmp (test, "GRAPPA")               ) (remote) ?   grappatest (rc) :   grappatest (lc);
		else if (!strcmp (test, "KTPoints")             ) (remote) ?      ktptest (rc) :      ktptest (lc);
		else if (!strcmp (test, "CompressedSensing")    ) (remote) ?       cstest (rc) :       cstest (lc);
		else if (!strcmp (test, "mxtest")               ) (remote) ?       mxtest (rc) :       mxtest (lc);
		else if (!strcmp (test, "nitest")               ) (remote) ?       nitest (rc) :       nitest (lc);
		else if (!strcmp (test, "fftwtest")             ) (remote) ?     fftwtest (rc) :     fftwtest (lc);
		else if (!strcmp (test, "dwttest")              ) (remote) ?      dwttest (rc) :      dwttest (lc);
		else if (!strcmp (test, "algotest")             ) (remote) ?     algotest (rc) :     algotest (lc);
		else if (!strcmp (test, "RelativeSensitivities")) (remote) ?     resetest (rc) :     resetest (lc);
		else if (!strcmp (test, "VDSpiral"))              (remote) ? vdspiraltest (rc) : vdspiraltest (lc);
		else                                              (remote) ? internaltest (rc) : internaltest (lc);

		if (remote) 
			delete rc;
		else 
			delete lc;
		
		return 0;

	} else
		
		return 1;	

}



bool init (int argc, char** argv) {
	
	cout << endl;
#ifdef VERSION
	cout << "jrrs "         << VERSION                                        << endl;
#else
	cout << "jrrs "         << endl;
#endif
	cout << "juelich remote reconstruction service "                          << endl;
#ifdef SVN_REVISION
	cout << "Test client "  << " [build " << SVN_REVISION << "]"              << endl;
#else
	cout << "Test client "  << endl;
#endif

    cout << "Copyright (C) 2010-2011"                                         << endl;
	cout << "Kaveh Vahedipour - k.vahedipour@fz-juelich.de"                   << endl;
	cout << "Juelich Research Centre"                                         << endl;
	cout << "Institute of Neuroscience and Medicine"                          << endl;
	cout << "Medical Imaging Physics"                                         << endl;
	cout << endl;

	Options *opt = new Options();
	
	opt->addUsage  ("Usage: testclt --name <name> [OPTIONS]");
	opt->addUsage  ("");
	opt->addUsage  (" -n, --name    Remote service name (for example: ReconService)");
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
	opt->setFlag   ("remote"  ,'r');
	opt->setFlag   ("pulses"  ,'p');

	opt->setOption ("name"   , 'n');
	opt->setOption ("test"   , 't');
	opt->setOption ("verbose", 'v');
	opt->setOption ("config" , 'c');
	opt->setOption ("data"   , 'd');
	opt->setOption ("base"   , 'b');
	
	opt->processCommandArgs(argc, argv);

	remote = opt->getFlag("remote");
	
	if ( !(opt->hasOptions()) || opt->getFlag("help") ) {
		
		opt->printUsage();
		delete opt;
		return 0;
		
	} 
	
	verbose = (opt->getValue("verbose" ) && atoi(opt->getValue("verbose" )) >= 0 && atoi(opt->getValue("verbose" )) <= 40) ?      opt->getValue("verbose" )  : (char*)"0";
	name    = (opt->getValue("name"  ) &&      opt->getValue("name"  )  != (char*)"")                                ?      opt->getValue("name"  )  : (char*)"ReconService" ;
	base    = (opt->getValue("base"  ) &&      opt->getValue("base"  )  != (char*)"")                                ?      opt->getValue("base"  )  : (char*)".";
	data    = (opt->getValue("data"  )     &&      opt->getValue("data"  )  != (char*)"")                                ?      opt->getValue("data"  )  : (char*)".";
	config  = (opt->getValue("config"  ) &&      opt->getValue("config"  )  != (char*)"")                                ?      opt->getValue("config"  )  : (char*)".";
	test    = (opt->getValue("test"  ) &&      opt->getValue("test"  )  != (char*)"")                                ?      opt->getValue("test"  )  : (char*)"DummyRecon";
	delete opt;

	return true;

}
