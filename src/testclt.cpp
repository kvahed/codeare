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

#include "options.h"
#include "ReconClient.hpp"

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
char*  debug;
char*  test;

bool init (int argc, char** argv);

int main (int argc, char** argv) {
	
	if (init (argc, argv)) {
		
		ReconClient client (name, debug);
		int         i = 0, j = 0, d = 8;
		
		Matrix<raw>    r (d, d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
		Matrix<double> h (d, d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
		Matrix<short>  p (d, d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);

		r.Random();
		h.Random();
		p.Random();

		std::cout << r << std::endl;
		std::cout << h << std::endl;
		std::cout << p << std::endl;
		
		client.SetRaw(r);
		client.SetPixel(p);
		client.SetHelper(h);
		client.SetKSpace(h);

		time_t seconds = time (NULL);
		char   uid[16];
		sprintf(uid,"%ld",seconds);
		
		client.ReadConfig("test.xml");
		client.SetAttribute("UID", uid);
		client.SetAttribute("Pi", 3.14156);
		client.SetAttribute("Dim", d);
		
		client.Process(test);

		client.GetRaw(r);
		client.GetPixel(p);
		client.GetHelper(h);
		client.GetKSpace(h);
		
		cout << "We're good" << endl;
		
		// Matrix<raw>    sensitivities;
		/*Matrix<raw>    data;
		Matrix<double> weights;
		Matrix<double> kspace;
		
		weights.read       ("share/nufft/data.h5", "weights");
		data.read          ("share/nufft/data.h5", "data");
		kspace.read        ("share/nufft/data.h5", "8_shot_piral");*/
		
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

    cout << "Copyright (C) 2010"                                              << endl;
	cout << "Kaveh Vahedipour - k.vahedipour@fz-juelich.de"                   << endl;
	cout << "Juelich Research Centre"                                         << endl;
	cout << "Institute of Neuroscience and Medicine"                          << endl;
	cout << "Medical Imaging Physics"                                         << endl;
	cout << endl;

	Options *opt = new Options();
	
	opt->addUsage  ("Usage: testclt --name <name> [OPTIONS]");
	opt->addUsage  ("");
	opt->addUsage  (" -n, --name   Remote service name (for example: ReconService)");
	opt->addUsage  (" -t, --test   Test case (default: DummyRecon. Just connectivity test)");
	opt->addUsage  (" -d, --debug  Debug level 0-40 (default: 0)");
	opt->addUsage  (" -b, --base   Base directory of approved files.");
	opt->addUsage  ("");
	opt->addUsage  (" -h, --help   Print this help screen"); 
	opt->addUsage  ("");
	
	opt->setFlag   ("help"  , 'h');

	opt->setOption ("name"  , 'n');
	opt->setOption ("test"  , 't');
	opt->setOption ("debug" , 'd');
	opt->setOption ("base"  , 'b');
	
	opt->processCommandArgs(argc, argv);
	
	if ( !(opt->hasOptions()) || opt->getFlag("help") ) {
		
		opt->printUsage();
		delete opt;
		return 0;
		
	} 
	
	debug = (opt->getValue("debug" ) && atoi(opt->getValue("debug" )) >= 0 && atoi(opt->getValue("debug" )) <= 40) ?      opt->getValue("debug" )  : (char*)"0";
	name  = (opt->getValue("name"  ) &&      opt->getValue("name"  )  != (char*)"")                                ?      opt->getValue("name"  )  : (char*)"ReconService" ;
	base  = (opt->getValue("base"  ) &&      opt->getValue("base"  )  != (char*)"")                                ?      opt->getValue("base"  )  : (char*)".";
	test  = (opt->getValue("test"  ) &&      opt->getValue("test"  )  != (char*)"")                                ?      opt->getValue("test"  )  : (char*)"DummyRecon";
	delete opt;

	return true;

}
