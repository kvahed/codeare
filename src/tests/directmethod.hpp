template <class T> bool 
dmtest (Connector<T>* rc) {

	Matrix<cxfl>  b1;  
	Matrix<cxfl>  tmxy;
	Matrix<float> tmz;  
	Matrix<float> r;   
	Matrix<float> b0;  
	
	Matrix<cxfl>  smxy;
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
	smxy = Matrix<cxfl>::Zeros  (r.Dim(1), 1);
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
	
	Matrix<cxfl>  rf;
	Matrix<cxfl>  mxy;
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

