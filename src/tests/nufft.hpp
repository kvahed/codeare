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
	IO::MXRead (weights, df, "weights");
	IO::MXRead (rawdata, df, "data");
	IO::MXRead  (kspace, df, "kspace");
#endif

	rc->SetMatrix    ("data",    rawdata);
	rc->SetMatrix    ("weights", weights);
	rc->SetMatrix    ("kspace",  kspace);
	
	rc->Process    (test);
	
	rc->GetMatrix    ("data", rawdata);

	rc->Finalise(test);

#ifdef HAVE_MAT_H	
	IO::MXDump   (rawdata, odf.c_str(), "image");
#endif
#ifdef HAVE_NIFTI1_IO_H
	IO::NIDump   (rawdata, "image.nii.gz");
#endif

	return true;
	
}
