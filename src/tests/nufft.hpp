template <class T> bool
nuffttest (Connector<T>* rc) {

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
