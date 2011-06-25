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
char*  data;
char*  config;
char*  verbose;
char*  test;

bool init (int argc, char** argv);

bool internaltest (ReconClient* rc); 
bool cgsensetest (ReconClient* rc);
bool nuffttest (ReconClient* rc);

int main (int argc, char** argv) {
	
	if (init (argc, argv)) {
		
		ReconClient client (name, verbose);

		if (strcmp (test, "NuFFT")   == 0)
			nuffttest (&client);
		else if (strcmp (test, "NuFFT_OMP")   == 0)
			nuffttest (&client);
		else if (strcmp (test, "CGSENSE") == 0)
			cgsensetest (&client);
		else
			internaltest (&client);
			
		return 0;

	} else
		
		return 1;

}

bool cgsensetest (ReconClient* rc) {

	Matrix<raw>    rawdata;
	Matrix<double> weights;
	Matrix<double> kspace;
	Matrix<raw>    sensitivities;
	
	std::string    cf  = std::string (base + std::string(config));

	std::string    df  = std::string (base + std::string(data));
	std::string    odf = std::string (base + std::string("/images.h5"));
	std::string    opf = std::string (base + std::string("/pulses.h5"));

	rc->ReadConfig  (cf.c_str());
	
	weights.read       (df, "weights");
	rawdata.read       (df, "data");
	kspace.read        (df, "kspace");
	sensitivities.read (df, "sensitivities");
	
	rc->SetRaw (rawdata);
	rc->SetRHelper (sensitivities);
	rc->SetHelper (weights);
	rc->SetKSpace (kspace);
	
	rc->Process ("CGSENSE");
	
	rc->GetRaw(rawdata);
	rc->GetRHelper(sensitivities);
	
	rawdata.dump(odf.c_str());
	sensitivities.dump(opf.c_str());
	
}

bool nuffttest (ReconClient* rc) {

	Matrix<raw>    data;
	Matrix<double> weights;
	Matrix<double> kspace;
	
	weights.read ("share/nufft/data.h5", "weights");
	data.read ("share/nufft/data.h5", "data");
	kspace.read ("share/nufft/data.h5", "8_shot_spiral");
	
	rc->ReadConfig ("share/nufft/config.xml");
	rc->SetRaw (data);
	rc->SetHelper (weights);
	rc->SetKSpace (kspace);
	
	rc->Process (test);
	
	rc->GetRaw(data);
	data.dump("share/nufft/recon.h5");
	
}

bool internaltest (ReconClient* rc) {

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
	
	rc->SetRaw(r);
	rc->SetRHelper(r);
	rc->SetPixel(p);
	rc->SetHelper(h);
	rc->SetKSpace(h);
	
	time_t seconds = time (NULL);
	char   uid[16];
	sprintf(uid,"%ld",seconds);
	
	rc->ReadConfig("test.xml");
	rc->SetAttribute("UID", uid);
	rc->SetAttribute("Pi", 3.14156);
	rc->SetAttribute("Dim", d);
	
	rc->Process(test);
	
	rc->GetRaw(r);
	rc->GetRHelper(r);
	rc->GetPixel(p);
	rc->GetHelper(h);
	rc->GetKSpace(h);
	
	cout << "We're good" << endl;
	
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
	opt->addUsage  ("");
	opt->addUsage  (" -h, --help    Print this help screen"); 
	opt->addUsage  ("");
	
	opt->setFlag   ("help"   , 'h');

	opt->setOption ("name"   , 'n');
	opt->setOption ("test"   , 't');
	opt->setOption ("verbose", 'v');
	opt->setOption ("config" , 'c');
	opt->setOption ("data"   , 'd');
	opt->setOption ("base"   , 'b');
	
	opt->processCommandArgs(argc, argv);
	
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
