/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipﬂour
 *                               Forschungszentrum J√ºlich, Germany
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
#include "ReconContext.hpp"
#include "FFT.hpp"

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

bool internaltest (ReconClient* rc); 
bool cgsensetest (ReconClient* rc);
bool nuffttest (ReconClient* rc);
bool sdmtest (ReconClient* rc);
bool mxtest (ReconClient* rc);
bool nitest (ReconClient* rc);
bool fftwtest (ReconClient* rc);
bool resetest (ReconClient* rc);
bool grappatest (ReconClient* rc);

int main (int argc, char** argv) {
	
	if (init (argc, argv)) {
		
		ReconClient client (name, verbose);

		if (strcmp (test, "NuFFT")   == 0)
			nuffttest (&client);
		else if (strcmp (test, "NuFFT_OMP")   == 0)
			nuffttest (&client);
		else if (strcmp (test, "CGSENSE") == 0)
			cgsensetest (&client);
		else if (strcmp (test, "GRAPPA") == 0)
			grappatest (&client);
		else if (strcmp (test, "SpatialDomain") == 0)
			sdmtest (&client);
		else if (strcmp (test, "mxtest") == 0)
			mxtest (&client);
		else if (strcmp (test, "nitest") == 0)
			nitest (&client);
		else if (strcmp (test, "fftwtest") == 0)
			fftwtest (&client);
		else if (strcmp (test, "RelativeSensitivities") == 0)
			resetest (&client);
		else
			internaltest (&client);
			
		return 0;

	} else
		
		return 1;

}

bool grappatest (ReconClient* rc) {

	Matrix<cplx> sig;
	Matrix<cplx> acs;
	
	std::string cf = std::string (base + std::string(config));
	std::string df = std::string (base + std::string(data));

#ifdef HAVE_MAT_H
	sig.MXRead  (df, "data", "");
	acs.MXRead  (df, "acs",  "");
#endif

	rc->ReadConfig (cf.c_str());
	
	rc->Init (test);
	rc->SetCplx  ("acs",  acs);
	rc->Prepare (test);
	rc->SetCplx  ("data", sig); // Measurement data
	rc->Process (test);
	rc->GetCplx  ("data", sig);
	rc->Finalise (test);
	
	return true;

}

bool cgsensetest (ReconClient* rc) {

	Matrix<cplx> rawdata;
	Matrix<double> weights;
	Matrix<double> kspace;
	Matrix<cplx> sens;
	
	std::string    cf  = std::string (base + std::string(config));
	std::string    df  = std::string (base + std::string(data));
	std::string    odf = std::string (base + std::string("/images.mat"));
	std::string    opf = std::string (base + std::string("/pulses.mat"));
	std::string    oif = std::string (base + std::string("/iters.mat"));

	weights.Read   (df, "weights");
	rawdata.Read   (df, "data");
	kspace.Read    (df, "kspace");
	sens.Read      (df, "sensitivities");

	if (remote) {
	
		rc->ReadConfig (cf.c_str());

		if (rc->Init (test) != OK) {
			printf ("Intialising failed ... bailing out!"); 
			return false;
		}
		
		// Outgoing -------------
		
		rc->SetCplx ("data", rawdata); // Measurement data
		rc->SetCplx ("sens", sens);    // Sensitivities
		rc->SetReal ("weights", weights); // Weights
		rc->SetReal ("kspace", kspace);  // K-space

		// ---------------------

		rc->Process    (test);
		
		// Incoming -------------

		rc->GetCplx     ("data", rawdata);  // Images
		if (pulses)
			rc->GetCplx ("sens", sens);     // Pulses (Excitation)
		rc->GetReal  ("weights", weights);  // CG residuals

		// ---------------------

		rc->Finalise   (test);

	} else {

		/*	RRServer::ReconContext* rx = new RRServer::ReconContext(test);

		rx->ReadConfig(cf.c_str());

		rx->SetCplx     (&rawdata);
		rawdata.Clear();

		rx->SetCplx (&sens);
		sens.Clear();

		rx->SetReal  (&weights);
		weights.Clear();

		rx->SetReal  (&kspace);
		kspace.Clear();

		rx->Init();
		rx->Process    ();
		
		rx->GetCplx     (&rawdata);
		rx->GetCplx (&sens);
		rx->GetReal  (&weights);
		*/			
	}
	
#ifdef HAVE_MAT_H	
	rawdata.MXDump   (odf.c_str(), "images");
	if (pulses)
		sens.MXDump  (opf.c_str(), "pulses");
	weights.MXDump   (oif.c_str(), "residuals");
#endif

	return true;
	
}

bool nuffttest (ReconClient* rc) {

	Matrix<cplx>   rawdata;
	Matrix<double> weights;
	Matrix<double> kspace;
	
	std::string    cf  = std::string (base + std::string(config));
	std::string    df  = std::string (base + std::string(data));
	std::string    odf = std::string (base + std::string("/images.mat"));

	rc->ReadConfig (cf.c_str());
	rc->Init(test);

	weights.MXRead   (df, "weights");
	rawdata.MXRead   (df, "data");
	kspace.MXRead    (df, "kspace");

	rc->SetCplx    ("data",    rawdata);
	rc->SetReal    ("weights", weights);
	rc->SetReal    ("kspace",  kspace);
	
	rc->Process    (test);
	
	rc->GetCplx    ("data", rawdata);

	rc->Finalise(test);

#ifdef HAVE_MAT_H	
	rawdata.MXDump   (odf.c_str(), "image");
#endif
#ifdef HAVE_NIFTI1_IO_H
	rawdata.NIDump   ("image.nii.gz");
#endif

	return true;
	
}


bool sdmtest (ReconClient* rc) {

	Matrix<cplx>   target;
	Matrix<cplx>   b1;
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

	rc->SetCplx    ("target", target);
	rc->SetCplx    ("b1",     b1);
	rc->SetReal    ("r",      r);
	rc->SetReal    ("k",      k);
	rc->SetPixel   ("b0",     b0);
	
	rc->Process    (test);
	
	rc->GetCplx    ("target", target);
	rc->GetCplx    ("b1",     b1);
	rc->GetReal    ("r",      r);

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

bool internaltest (ReconClient* rc) {

	int            i = 0, j = 0, d = 5;
	
	Matrix<cplx>   r (d, d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
	Matrix<double> h (d, d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
	Matrix<short>  p (d, d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
	
	r.Random(); 
	h.Random();
	p.Random();

	Matrix<size_t> m (3,2);
	m (0,0) = 1;
	m (0,1) = 3;
	m (1,0) = 1;
	m (1,1) = 4;
	m (2,0) = 1;
	m (2,1) = 5;

	std::cout << m << std::endl;
	
	Matrix<size_t> mg = Matrix<size_t>::MeshGrid(m);
	mg.MXDump("mg.mat", "mg");
	
	std::cout << r << std::endl;
	std::cout << h << std::endl;
	std::cout << p << std::endl;
	
	Matrix<std::complex<double> >  f = (Matrix<std::complex<double> >) r;

	rc->ReadConfig("test.xml");
	rc->Init(test);

	rc->SetCplx ("r", r);
	rc->SetPixel("p", p);
	rc->SetReal ("h", h);
	
	time_t seconds = time (NULL);
	char   uid[16];
	sprintf(uid,"%ld",seconds);
	
	rc->SetAttribute("UID", uid);
	rc->SetAttribute("Pi", 3.14156);
	rc->SetAttribute("Dim", d);
	
	rc->Process(test);
	
	rc->GetCplx ("r", r);
	rc->GetPixel("p", p);
	rc->GetReal ("h", h);

	rc->Finalise (test);
	
	cout << "We're good" << endl;

	return true;
	
}

bool fftwtest (ReconClient* rc) {

	std::string in  = std::string (base + std::string ("/infft.h5"));
	std::string out = std::string (base + std::string ("/outfft.h5"));
	
	Matrix<cplx> m;

	m.Read (in, "img");
	m = FFT::Forward(m);
	m = m.FFTShift();
	m = m.IFFTShift();
	m = FFT::Backward(m);
	m.Dump (out, "img");

	return true;

}

bool resetest (ReconClient* rc) {

	// OUT:
	Matrix<cplx>   meas; // measurement
	Matrix<cplx>   mask; // measurement

	// IN: 
	Matrix<cplx>   txm;  // Transmit maps
	Matrix<cplx>   rxm;  // Receive maps

	Matrix<double> b0;   // B0 map
	Matrix<double> snro; // SNR optimal image
	Matrix<double> bet;  // GRE image mask

	Matrix<short> bets;  // BET mask
	// Read configuration file and initialise backend ------------

	std::string    cf  = std::string (base + std::string (config));
	rc->ReadConfig (cf.c_str());
	rc->DumpConfig ("con.xml");
	
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
	mask.RAWRead (maf, std::string("VB15"));

	rc->SetCplx ("meas", meas);
	rc->SetCplx ("mask", mask);
	// -----------------------------------------------------------

	// Process data on backend -----------------------------------

	rc->Process(test);
	// -----------------------------------------------------------
	
	// Get back reconstructed data from backend ------------------

	rc->GetCplx ("txm",  txm);
	rc->GetCplx ("rxm",  rxm);
	rc->GetCplx ("mask", mask);
	rc->GetReal ("snro", snro);
	rc->GetReal ("b0",   b0);
	rc->GetPixel("bets", bets);

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

bool mxtest (ReconClient* rc) {

#ifdef HAVE_MAT_H

	Matrix<double> in (3,5);
	in.Random ();

	std::cout << in << std::endl << std::endl;

	in.MXDump(std::string("test.mat"), std::string("imat"), std::string(""));

	Matrix<double> out;
	out.MXRead(std::string("test.mat"), std::string("imat"), std::string(""));

	std::cout << out << std::endl << std::endl;

	Matrix<cplx> r1 (4,8);
	r1.Random ();
	r1.MXDump (std::string("rtest.mat"), std::string("rmat"), std::string(""));
	std::cout << r1 << std::endl << std::endl;
	
	Matrix<cplx> r2;
	r2.MXRead (std::string("rtest.mat"), std::string("rmat"), std::string(""));
	std::cout << r2 << std::endl << std::endl;

#else

	std::cout << "MATLAB root not set during configuration (--with-matlabroot).\n Test skipped." << std::endl;

#endif

	return true;

}

bool nitest (ReconClient* rc) {

#ifdef HAVE_NIFTI1_IO_H
	
	std::string    cf  = std::string (base + std::string(config));
	std::string    df  = std::string (base + std::string(data));
	std::string    mat = std::string (base + std::string("betted.mat"));
	std::string    nii = std::string (base + std::string("betted2.nii.gz"));
	
	Matrix<double> d;
	d.NIRead (df);
	d.MXDump (mat, std::string("betted"), std::string(""));
	d.NIDump (nii);

	Matrix<cplx> slp = Matrix<cplx>::Phantom3D(196); 
	slp.NIDump (nii);
	slp.MXDump (mat, std::string("betted"), std::string(""));

#else

	std::cout << "No nifti support compiled in. Bailing out." << std::endl;

#endif

	return true;

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
