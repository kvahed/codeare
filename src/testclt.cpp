#include "options.h"
#include "Matrix.h"
#include "ReconClient.h"
#ifndef __WIN32__
#include "config.h"
#endif

#ifndef SVN_REVISION
	#define SVN_REVISION "unkown"
#endif

using namespace std;

char*  name;
char*  base;
char*  debug;
int    test;

bool init (int argc, char** argv);


int main (int argc, char** argv) {
	
	if (init (argc, argv)) {

		ReconClient client (name, debug);
		
		int         i, j, d;
		
		d = 2;
		Matrix< complex<float> > data (d, d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
		
		for (i = 0; i < d; i++)
			for (j = 0; j < d; j++)
				data.at(i,j) = complex<float> ((float) i, (float) j);
		
		client.SetRaw(data);
		client.Requestprocess_data((RRSModule::method) test);
		client.GetRaw(data); 
		
		for (i = 0; i < d; i++)
			for (j = 0; j < d; j++)
				cout << data.at(i,j) << endl;

		Matrix< short > pdata (d, d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
		
		for (i = 0; i < d; i++)
			for (j = 0; j < d; j++)
				pdata.at(i,j) = (i+1)*(j+1);
		
		client.SetPixel(pdata);
		client.Requestprocess_data((RRSModule::method) test);
		client.GetPixel(pdata); 
		
		for (i = 0; i < d; i++)
			for (j = 0; j < d; j++)
				cout << pdata.at(i,j) << endl;

		return 0;

	} else
		
		return 1;

}


bool init (int argc, char** argv) {

	cout << endl;
	cout << "jrrs "         << VERSION                                        << endl;
	cout << "juelich remote reconstruction service "                          << endl;
	cout << "Test client "  << " [build " << SVN_REVISION << "]"              << endl;
    cout << "Copyright (C) 2010	Kaveh Vahedipour<k.vahedipour@fz-jeulich.de>" << endl;
	cout << "Kaveh Vahedipour -  k.vahedipour@fz-juelich.de"                  << endl;
	cout << "Juelich Research Centre"                                         << endl;
	cout << "Institute of Neuroscience and Medicine"                          << endl;
	cout << "Medical Imaging Physics"                                         << endl;
	cout << endl;

	Options *opt = new Options();
	
	opt->addUsage  ("Usage: testclt [OPTIONS]");
	opt->addUsage  ("");
	opt->addUsage  (" -n, --name   Remote service name (default: ReconService)");
	opt->addUsage  (" -t, --test   Test case (default: 0, no recon. Just connectivity test)");
	opt->addUsage  (" -d, --debug  Debug level 0-40 (default: 0)");
	opt->addUsage  (" -b, --base   Base directory of approved files.");
	opt->addUsage  ("");
	opt->addUsage  (" -h, --help   Print this help screen"); 
	opt->addUsage  ("");
	opt->addUsage  ("If no options are given the following options are assumed:");
	opt->addUsage  ("testclt -n ReconService -t 0 -d 0 -b ." );
	opt->addUsage  ("");
	
	
	opt->setFlag   ("help"  , 'h');

	opt->setOption ("name"  , 'n');
	opt->setOption ("test"  , 't');
	opt->setOption ("debug" , 'd');
	opt->setOption ("base"  , 'b');
	
	opt->processCommandArgs(argc, argv);
	
	if ( opt->getFlag("help") ) {
		
		opt->printUsage();
		delete opt;
		return 0;
		
	} 
	
	test  = (opt->getValue("test"  ) && atoi(opt->getValue("test"  )) >= 0)                                        ? atoi(opt->getValue("test"  )) : 0 ;
	debug = (opt->getValue("debug" ) && atoi(opt->getValue("debug" )) >= 0 && atoi(opt->getValue("debug" )) <= 40) ?      opt->getValue("debug" )  : (char*)"0";
	name  = (opt->getValue("name"  ) &&      opt->getValue("name"  )  != (char*)"")                                ?      opt->getValue("name"  )  : (char*)"ReconService" ;
	base  = (opt->getValue("base"  ) &&      opt->getValue("base"  )  != (char*)"")                                ?      opt->getValue("base"  )  : (char*)".";

	delete opt;

	return true;

}
