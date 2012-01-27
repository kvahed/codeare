template <class T> bool 
internaltest (Connector<T>* rc) {

	int            i = 0, j = 0, d = 5;
	
	Matrix<cxfl>   r (d, d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
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

