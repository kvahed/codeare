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

	IO::Read  (target, df, "target");
	IO::Read      (b1, df, "b1");
	IO::Read      (b0, df, "b0");
	IO::Read       (k, df, "k");
	IO::Read       (r, df, "r");

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

	IO::MXDump (target, mf, "pattern");
	IO::MXDump     (b1, mf, "ptx");
	IO::MXDump      (r, mf, "nrmse");

	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n",fname.c_str());
		return false;
	}
#endif
	return true;

}

