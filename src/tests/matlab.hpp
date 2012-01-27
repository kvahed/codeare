template <class T> bool
mxtest (Connector<T>* rc) {

#ifdef HAVE_MAT_H

	Matrix<double> in (3,5);
	in.Random ();

	std::cout << in << std::endl << std::endl;

	in.MXDump(std::string("test.mat"), std::string("imat"), std::string(""));

	Matrix<double> out;
	out.MXRead(std::string("test.mat"), std::string("imat"), std::string(""));

	std::cout << out << std::endl << std::endl;

	Matrix<cxfl> r1 (4,8);
	r1.Random ();
	r1.MXDump (std::string("rtest.mat"), std::string("rmat"), std::string(""));
	std::cout << r1 << std::endl << std::endl;
	
	Matrix<cxfl> r2;
	r2.MXRead (std::string("rtest.mat"), std::string("rmat"), std::string(""));
	std::cout << r2 << std::endl << std::endl;

#else

	std::cout << "MATLAB root not set during configuration (--with-matlabroot).\n Test skipped." << std::endl;

#endif

	return true;

}

