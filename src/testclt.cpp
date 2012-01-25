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

#include "options.h"
#include "RemoteConnector.hpp"
#include "LocalConnector.hpp"
#include "modules/FFT.hpp"
#include "MatrixOperations.hpp"
//#include "tests/cgsense.hpp"

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
bool   remote;
bool   pulses;


bool init (int argc, char** argv);

template <class T> bool 
internaltest (Connector<T>* rc); 

template <class T> bool 
cgsensetest (Connector<T>* rc); 

/*
void selector ()
{
   float result = pt2Func(a, b);    // call using function pointer

   cout << "Switch replaced by function pointer: 2-5=";  // display result
   cout << result << endl;
}
*/

int main (int argc, char** argv) {
	
	if (init (argc, argv)) {

		
		Connector<RemoteConnector>* rcon; 
		Connector<LocalConnector>*  lcon;

		if (remote) 
			rcon = new Connector<RemoteConnector> (name, verbose);
		else 
			lcon = new Connector<LocalConnector> (name, verbose);
		
		if (strcmp      (test, "CGSENSE")   == 0)
			(remote) ? cgsensetest (rcon) : cgsensetest (lcon);
		/*
		  else if (strcmp      (test, "NuFFT")   == 0)
		  nuffttest (conn);
		  else if (strcmp (test, "NuFFT_OMP")   == 0)
		  nuffttest (conn);
		  else if (strcmp (test, "GRAPPA") == 0)
		  grappatest (conn);
		  else if (strcmp (test, "KTPoints") == 0)
		  ktptest (conn);
		  else if (strcmp (test, "CompressedSensing") == 0)
		  cstest (conn);
		  else if (strcmp (test, "mxtest") == 0)
		  mxtest (conn);
		  else if (strcmp (test, "nitest") == 0)
		  nitest (conn);
		  else if (strcmp (test, "fftwtest") == 0)
		  fftwtest (conn);
		  else if (strcmp (test, "RelativeSensitivities") == 0)
		  resetest (conn);
		  else if (strcmp (test, "DirectMethod") == 0)
		  dmtest (conn);
		  else
		  (remote) ? internaltest (rc) : internaltest (lc);
		*/
		
		return 0;

		if (remote) 
			delete rcon;
		else 
			delete lcon;
		
	} else
		
		return 1;	

}

template <class T> bool 
cgsensetest (RRClient::Connector<T>* rc) {

	// Incoming
	Matrix<cxfl>   rawdata;
	Matrix<double> weights;
	Matrix<double> kspace;
	Matrix<cxfl>   sens;
	
	// Outgoing
	Matrix<double> nrmse;
	Matrix<cxfl>   image;
	Matrix<cxfl>   signals;
	
	std::string    cf  = std::string (base + std::string(config));
	std::string    df  = std::string (base + std::string(data));
	std::string    odf = std::string (base + std::string("/images.mat"));
	
	weights.Read   (df, "weights");
	rawdata.Read   (df, "data");
	kspace.Read    (df, "kspace");
	sens.Read      (df, "sensitivities");
	
	rc->ReadConfig (cf.c_str());
	
	if (rc->Init (test) != OK) {
		printf ("Intialising failed ... bailing out!\n"); 
		return false;
	}
	
	rc->Attribute ("pulses", (int*)&pulses);
	
	// Outgoing -------------
	
	rc->SetMatrix (   "data", rawdata); // Measurement data
	rc->SetMatrix (   "sens", sens);    // Sensitivities
	rc->SetMatrix ("weights", weights); // Weights
	rc->SetMatrix ( "kspace", kspace);  // K-space
	
	// ---------------------
	
	rc->Process    (test);
	
	// Incoming -------------
	
	rc->GetMatrix     (  "image", image);  // Images
	if (pulses)
		rc->GetMatrix ("signals", signals);     // Pulses (Excitation)
	rc->GetMatrix     (  "nrmse", nrmse);  // CG residuals
	
	// ---------------------
	
	rc->Finalise   (test);
	
#ifdef HAVE_MAT_H	
	MATFile* mf = matOpen (odf.c_str(), "w");
	
	if (mf == NULL) {
		printf ("Error creating file %s\n", odf.c_str());
		return false;
	}
	
	image.MXDump       (mf, "image");
	if (pulses)
		signals.MXDump (mf, "signals");
	if (nrmse.Size() > 1)
		nrmse.MXDump   (mf, "nrmse");
	
	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", odf.c_str());
		return false;
	}
#endif

	return true;
	
}

/*
bool nuffttest (Connector* rc) {

	Matrix<cxfl>   rawdata;
	Matrix<double> weights;
	Matrix<double> kspace;
	
	std::string    cf  = std::string (base + std::string(config));
	std::string    df  = std::string (base + std::string(data));
	std::string    odf = std::string (base + std::string("/images.mat"));

	rc->ReadConfig (cf.c_str());
	rc->Init(test);

#ifdef HAVE_MAT_H	
	weights.MXRead (df, "weights");
	rawdata.MXRead (df, "data");
	kspace.MXRead  (df, "kspace");
#endif

	rc->SetMatrix    ("data",    rawdata);
	rc->SetMatrix    ("weights", weights);
	rc->SetMatrix    ("kspace",  kspace);
	
	rc->Process    (test);
	
	rc->GetMatrix    ("data", rawdata);

	rc->Finalise(test);

#ifdef HAVE_MAT_H	
	rawdata.MXDump   (odf.c_str(), "image");
#endif
#ifdef HAVE_NIFTI1_IO_H
	rawdata.NIDump   ("image.nii.gz");
#endif

	return true;
	
}


bool ktptest (Connector* rc) {

	Matrix<cxfl>   target;
	Matrix<cxfl>   b1;
	Matrix<double> r;
	Matrix<double> k;
	Matrix<short>  b0;
	
	std::string    cf  = std::string (base + std::string(config));
	std::string    df  = std::string (base + std::string(data));
	std::string    odf = std::string (base + std::string("sdmout.mat"));
	std::string    pdf = std::string (base + std::string("result.h5"));
	std::string    rdf = std::string (base + std::string("residuals.h5"));

	rc->ReadConfig (cf.c_str());
	rc->Init(test);

	target.Read  (df, "target");
	b1.Read      (df, "b1");
	b0.Read      (df, "b0");
	k.Read       (df, "k");
	r.Read       (df, "r");

	rc->SetMatrix    ("target", target);
	rc->SetMatrix    ("b1",     b1);
	rc->SetMatrix    ("r",      r);
	rc->SetMatrix    ("k",      k);
	rc->SetMatrix   ("b0",     b0);
	
	rc->Process    (test);
	
	rc->GetMatrix    ("target", target);
	rc->GetMatrix    ("b1",     b1);
	rc->GetMatrix    ("r",      r);

	rc->Finalise(test);
	
	std::string fname = std::string (base + std::string ("sdout.mat"));
	
#ifdef HAVE_MAT_H	
	MATFile* mf = matOpen (fname.c_str(), "w");

	if (mf == NULL) {
		printf ("Error creating file %s\n", fname.c_str());
		return false;
	}

	target.MXDump (mf, "pattern", "");
	b1.MXDump     (mf, "ptx", "");
	r.MXDump      (mf, "nrmse", "");

	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n",fname.c_str());
		return false;
	}
#endif
	return true;

}

bool fftwtest (Connector* rc) {

	std::string in  = std::string (base + std::string ("/infft.h5"));
	std::string out = std::string (base + std::string ("/outfft.h5"));
	
	Matrix<cxfl> m;

	m.Read (in, "img");
	m = FFT::Forward(m);
	m = FFT::Backward(m);
	m.Dump (out, "img");

	return true;

}

bool resetest (Connector* rc) {

	// OUT:
	Matrix<cxfl>   meas; // measurement
	Matrix<cxfl>   mask; // measurement

	// IN: 
	Matrix<cxfl>   txm;  // Transmit maps
	Matrix<cxfl>   rxm;  // Receive maps

	Matrix<double> b0;   // B0 map
	Matrix<double> snro; // SNR optimal image
	Matrix<double> bet;  // GRE image mask

	Matrix<short> bets;  // BET mask
	// Read configuration file and initialise backend ------------

	std::string    cf  = std::string (base + std::string (config));
	rc->ReadConfig (cf.c_str());

	int use_bet = 0;
	rc->Attribute ("use_bet", &use_bet);

	stringstream ss;
	string mef, maf;

	ss << base << rc->Attribute("meas");
	mef = ss.str();
	ss.str("");
	ss << base << rc->Attribute("mask");
	maf = ss.str();

	rc->Init(test);
	// -----------------------------------------------------------

	// Read binary data and transmit to backend ------------------ 

	//sprintf ("--- %s ---\n", mef.c_str());

	meas.RAWRead (mef, std::string("VB15"));
	if (use_bet==1)
		mask.RAWRead (maf, std::string("VB15"));

	rc->SetMatrix ("meas", meas);
	rc->SetMatrix ("mask", mask);
	// -----------------------------------------------------------

	// Process data on backend -----------------------------------

	rc->Process(test);
	// -----------------------------------------------------------
	
	// Get back reconstructed data from backend ------------------

	rc->GetMatrix ("txm",  txm);
	rc->GetMatrix ("rxm",  rxm);
	rc->GetMatrix ("mask", mask);
	rc->GetMatrix ("snro", snro);
	rc->GetMatrix ("b0",   b0);
	rc->GetMatrix("bets", bets);

	// -----------------------------------------------------------

	// Clear RAM and hangup --------------------------------------

	rc->Finalise(test);
	// -----------------------------------------------------------

	// Write data to a single matlab disk ------------------------
	
	std::string fname = std::string (base + std::string ("maps.mat"));
	
#ifdef HAVE_MAT_H	
	MATFile* mf = matOpen (fname.c_str(), "w");

	if (mf == NULL) {
		printf ("Error creating file %s\n", fname.c_str());
		return false;
	}

	txm.MXDump  (mf,  "txm", "");
	rxm.MXDump  (mf,  "rxm", "");
	snro.MXDump (mf, "snro", "");
	b0.MXDump   (mf,   "b0", "");
	bets.MXDump (mf, "bets", "");
	mask.MXDump (mf, "mask", "");

	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n",fname.c_str());
		return false;
	}
#endif
	// -----------------------------------------------------------

	return true;

}

bool mxtest (Connector* rc) {

#ifdef HAVE_MAT_H

	Matrix<double> in (3,5);
	in.Random ();

	std::cout << in << std::endl << std::endl;

	in.MXDump(std::string("test.mat"), std::string("imat"), std::string(""));

	Matrix<double> out;
	out.MXRead(std::string("test.mat"), std::string("imat"), std::string(""));

	std::cout << out << std::endl << std::endl;

	Matrix<cxfl> r1 (4,8);
	r1.Random ();
	r1.MXDump (std::string("rtest.mat"), std::string("rmat"), std::string(""));
	std::cout << r1 << std::endl << std::endl;
	
	Matrix<cxfl> r2;
	r2.MXRead (std::string("rtest.mat"), std::string("rmat"), std::string(""));
	std::cout << r2 << std::endl << std::endl;

#else

	std::cout << "MATLAB root not set during configuration (--with-matlabroot).\n Test skipped." << std::endl;

#endif

	return true;

}

bool nitest (Connector* rc) {

#ifdef HAVE_NIFTI1_IO_H
	
	std::string    cf  = std::string (base + std::string(config));
	std::string    df  = std::string (base + std::string(data));
	std::string    mat = std::string (base + std::string("betted.mat"));
	std::string    nii = std::string (base + std::string("betted2.nii.gz"));
	
	Matrix<double> d;
	d.NIRead (df);
#ifdef HAVE_MAT_H	
	d.MXDump (mat, std::string("betted"), std::string(""));
#endif
	d.NIDump (nii);

	Matrix<cxfl> slp = Matrix<cxfl>::Phantom3D(196); 
	slp.NIDump (nii);
#ifdef HAVE_MAT_H	
	slp.MXDump (mat, std::string("betted"), std::string(""));
#endif

#else

	std::cout << "No nifti support compiled in. Bailing out." << std::endl;

#endif

	return true;

}
*/

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
