template <class T> bool
ktptest (Connector<T>* rc) {

	Matrix<cxfl>   target;
	Matrix<cxfl>   b1;
	Matrix<cxfl>   final;
	Matrix<float>  r;
	Matrix<float>  k;
	Matrix<float>  b0;
	
	std::string    cf  = std::string (base + std::string(config));
	std::string    df  = std::string (base + std::string(data));
	std::string    odf = std::string (base + std::string("sdmout.mat"));

	rc->ReadConfig (cf.c_str());
	rc->Init(test);

#ifdef HAVE_MAT_H
    #define FORMAT MATLAB
#else
	#define FORMAT HDF5
#endif

	Read  (target, df, "target", "", FORMAT);
	Read  (b1    , df, "b1"    , "", FORMAT);
	Read  (b0    , df, "b0"    , "", FORMAT);
	Read  (k     , df, "k"     , "", FORMAT);
	Read  (r     , df, "r"     , "", FORMAT);

	rc->SetMatrix ("target", target);
	rc->SetMatrix ("b1",     b1);
	rc->SetMatrix ("r",      r);
	rc->SetMatrix ("k",      k);
	rc->SetMatrix ("b0",     b0);
	
	rc->Process    (test);

	Matrix<cxfl>   rf;
	Matrix<float>  grad;
	Matrix<float>  nrmse;

	rc->GetMatrix ("nrmse",  nrmse);
	rc->GetMatrix ("rf",     rf);
	rc->GetMatrix ("grad",   grad);
	rc->GetMatrix ("final", final);

	rc->Finalise(test);
	
	std::string fname = std::string (base + std::string ("sdout.mat"));
	
#ifdef HAVE_MAT_H	
	MATFile* mf = matOpen (fname.c_str(), "w");

	if (mf == NULL) {
		printf ("Error creating file %s\n", fname.c_str());
		return false;
	}

	MXDump (nrmse, mf, "nrmse");
	MXDump ( grad, mf,  "grad");
	MXDump (   rf, mf,    "rf");
	MXDump (final, mf, "final");

	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n",fname.c_str());
		return false;
	}
#endif

	return true;

}

