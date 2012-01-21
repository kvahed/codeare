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
#include "ReconClient.hpp"
//#include "ReconContext.hpp"
#include "modules/FFT.hpp"
#include "MatrixOperations.hpp"

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
bool cgsensetest  (ReconClient* rc);
bool nuffttest    (ReconClient* rc);
bool ktptest      (ReconClient* rc);
bool mxtest       (ReconClient* rc);
bool nitest       (ReconClient* rc);
bool fftwtest     (ReconClient* rc);
bool resetest     (ReconClient* rc);
bool grappatest   (ReconClient* rc);
bool cstest       (ReconClient* rc);
bool dmtest       (ReconClient* rc);

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
		else if (strcmp (test, "KTPoints") == 0)
			ktptest (&client);
		else if (strcmp (test, "CompressedSensing") == 0)
			cstest (&client);
		else if (strcmp (test, "mxtest") == 0)
			mxtest (&client);
		else if (strcmp (test, "nitest") == 0)
			nitest (&client);
		else if (strcmp (test, "fftwtest") == 0)
			fftwtest (&client);
		else if (strcmp (test, "RelativeSensitivities") == 0)
			resetest (&client);
		else if (strcmp (test, "DirectMethod") == 0)
			dmtest (&client);
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
	rc->SetMatrix  ("acs",  acs);
	rc->Prepare (test);
	rc->SetMatrix  ("data", sig); // Measurement data
	rc->Process (test);
	rc->GetMatrix  ("data", sig);
	rc->Finalise (test);
	
	return true;

}

bool cgsensetest (ReconClient* rc) {

	// Incoming
	Matrix<cplx>   rawdata;
	Matrix<double> weights;
	Matrix<double> kspace;
	Matrix<cplx>   sens;

	// Outgoing
	Matrix<double> nrmse;
	Matrix<cplx>   image;
	Matrix<cplx>   signals;
	
	std::string    cf  = std::string (base + std::string(config));
	std::string    df  = std::string (base + std::string(data));
	std::string    odf = std::string (base + std::string("/images.mat"));

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
		
	} else {
		/*
		RRServer::ReconContext* rx = new RRServer::ReconContext(test);

		rx->ReadConfig(cf.c_str());
		rx->Init();

		rx->SetMatrix ("data",     rawdata);
		rx->SetMatrix ("sens",     sens);
		rx->SetMatrix ("weights",  weights);
		rx->SetMatrix ("kspace",   kspace);

		rx->Process ();
		
		rx->GetMatrix ("data",     rawdata); // Images
		if (pulses)
			rx->GetMatrix ("sens", sens);    // Pulses (Excitation)
		rx->GetMatrix ("weights",  weights); // CG residuals

		rx->Finalise ();
		*/
	}
	
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

bool cstest (ReconClient* rc) {

	Matrix<cplx> indata;
	Matrix<cplx> im_dc;
	Matrix<double> mask;
	Matrix<double> pdf;
	
	std::string cf  = std::string (base + std::string(config));
	std::string df  = std::string (base + std::string(data));
	std::string odf = std::string (base + std::string("/csout.mat"));

#ifdef HAVE_MAT_H	
	indata.MXRead (df, "data");
	pdf.MXRead (df, "pdf");
	mask.MXRead (df, "mask");
#endif

	rc->ReadConfig (cf.c_str());
	
	if (rc->Init (test) != OK) {
		printf ("Intialising failed ... bailing out!"); 
		return false;
	}
	
	// Outgoing -------------
	
	rc->SetMatrix  ("data", indata); // Measurement data
	rc->SetMatrix  ("pdf",  pdf);  // Sensitivities
	rc->SetMatrix  ("mask", mask); // Weights
	
	// ---------------------
	
	rc->Process (test);
	
	// Incoming -------------
	
	rc->GetMatrix ("data", indata);  // Images
	rc->GetMatrix ("im_dc", im_dc);  // Images
	
	// ---------------------
	
	rc->Finalise   (test);
	
#ifdef HAVE_MAT_H	
	MATFile* mf = matOpen (odf.c_str(), "w");

	if (mf == NULL) {
		printf ("Error creating file %s\n", odf.c_str());
		return false;
	}

	indata.MXDump (mf, "img");
	im_dc.MXDump (mf, "wvt");

	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", odf.c_str());
		return false;
	}
#endif
	
	return true;
	
}

bool dmtest (ReconClient* rc) {

	Matrix<cplx>  b1;  
	Matrix<cplx>  tmxy;
	Matrix<float> tmz;  
	Matrix<float> r;   
	Matrix<float> b0;  
	
	Matrix<cplx>  smxy;
	Matrix<float> smz;
	Matrix<float> roi;

	Matrix<float> g;   
	Matrix<float> j;
	
	std::string cf  = std::string (base + std::string(config));
	std::string odf = std::string (base + std::string("/simout.mat"));

	rc->ReadConfig (cf.c_str());
	
	std::string gf = std::string(base + std::string(rc->Attribute("gf"))); // gradient trajectories
	std::string pf = std::string(base + std::string(rc->Attribute("pf"))); // patterns
	std::string mf = std::string(base + std::string(rc->Attribute("mf"))); // maps

	// Gradients
#ifdef HAVE_MAT_H
	g.MXRead    (gf, rc->Attribute("g"));
	j.MXRead    (gf, rc->Attribute("j"));

	// Target excitation, ROI, sample
	r.MXRead    (pf, "r");
	tmxy.MXRead (pf, rc->Attribute("p"));
	tmz  = Matrix<float>::Zeros (r.Dim(1), 1);
	smxy = Matrix<cplx>::Zeros  (r.Dim(1), 1);
	smz.MXRead (pf, rc->Attribute("s"));
	roi.MXRead (pf, rc->Attribute("roi"));

	// Maps
	b1.MXRead (mf, rc->Attribute("b1"));
	b0.MXRead (mf, rc->Attribute("b0"));
#endif	
	if (rc->Init (test) != OK) {
		printf ("Intialising failed ... bailing out!"); 
		return false;
	}

	// Outgoing -------------
	
	rc->SetMatrix  (  "b1", b1  );
	rc->SetMatrix (   "g", g  );
	rc->SetMatrix (   "r", r   );
	rc->SetMatrix (  "b0", b0  );
	rc->SetMatrix  ("tmxy", tmxy);
	rc->SetMatrix ( "tmz", tmz );
	rc->SetMatrix  ("smxy", smxy);
	rc->SetMatrix ( "smz", smz );
	rc->SetMatrix (   "j", j   );
	rc->SetMatrix ( "roi", roi );
	// ---------------------
	
	rc->Process (test);
	
	// Incoming -------------
	
	Matrix<cplx>  rf;
	Matrix<cplx>  mxy;
	Matrix<float> mz;   

	rc->GetMatrix  ( "mxy", mxy);
	rc->GetMatrix (  "mz", mz);	
	rc->GetMatrix  ("tmxy", tmxy);
	rc->GetMatrix ( "tmz", tmz);	
	rc->GetMatrix  (  "rf", rf);

	// ---------------------
	
	rc->Finalise   (test);
	
#ifdef HAVE_MAT_H	
	MATFile* od = matOpen (odf.c_str(), "w");

	if (od == NULL) {
		printf ("Error creating file %s\n", odf.c_str());
		return false;
	}

	mxy.MXDump  (od, "mxy");
	mz.MXDump   (od, "mz");
	tmxy.MXDump (od, "tmxy");
	tmz.MXDump  (od, "tmz");
	rf.MXDump   (od, "rf");

	if (matClose(od) != 0) {
		printf ("Error closing file %s\n", odf.c_str());
		return false;
	}
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


bool ktptest (ReconClient* rc) {

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
#ifdef HAVE_MAT_H	
	mg.MXDump("mg.mat", "mg");
#endif
	
	std::cout << r << std::endl;
	std::cout << h << std::endl;
	std::cout << p << std::endl;
	
	Matrix<std::complex<double> >  f = (Matrix<std::complex<double> >) r;

	rc->ReadConfig("test.xml");
	rc->Init(test);

	rc->SetMatrix ("r", r);
	rc->SetMatrix ("p", p);
	rc->SetMatrix ("h", h);
	
	time_t seconds = time (NULL);
	char   uid[16];
	sprintf(uid,"%ld",seconds);
	
	rc->SetAttribute("UID", uid);
	rc->SetAttribute("Pi", 3.14156);
	rc->SetAttribute("Dim", d);
	
	rc->Process(test);
	
	rc->GetMatrix ("r", r);
	rc->GetMatrix("p", p);
	rc->GetMatrix ("h", h);

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
#ifdef HAVE_MAT_H	
	d.MXDump (mat, std::string("betted"), std::string(""));
#endif
	d.NIDump (nii);

	Matrix<cplx> slp = Matrix<cplx>::Phantom3D(196); 
	slp.NIDump (nii);
#ifdef HAVE_MAT_H	
	slp.MXDump (mat, std::string("betted"), std::string(""));
#endif

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
