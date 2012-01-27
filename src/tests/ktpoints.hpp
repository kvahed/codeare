template <class T> bool
ktptest (Connector<T>* rc) {

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

