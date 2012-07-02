template <class T> bool
nuffttest (Connector<T>* rc) {

	Matrix<cxdb>   rawdata;
	Matrix<double> weights;
	Matrix<double> kspace;
	Matrix<cxdb>   img;
	
	std::string    cf  = std::string (base + std::string(config));
	std::string    df  = std::string (base + std::string(data));
	std::string    odf = std::string (base + std::string("/images.mat"));

	rc->ReadConfig (cf.c_str());
	rc->Init(test);

#ifdef HAVE_MAT_H	
	if (!(MXRead (weights, df, "weights") && 
		  MXRead (rawdata, df, "data") && 
		  MXRead  (kspace, df, "kspace")))
		return false;
#endif
	
	rc->SetMatrix    ("data",    rawdata);
	rc->SetMatrix    ("weights", weights);
	rc->SetMatrix    ("kspace",  kspace);
	
	rc->Prepare      (test);
	rc->Process      (test);
	
	rc->GetMatrix    ("img", img);

	rc->Finalise(test);

#ifdef HAVE_MAT_H	
	MXDump   (img, odf.c_str(), "image");
#endif

	/*
#ifdef HAVE_NIFTI1_IO_H
	NIDump   (img, "image.nii.gz");
#endif
	*/
	return true;
	
}
