template <class T> bool
mxtest (Connector<T>* rc) {

#ifdef HAVE_MAT_H

	Matrix<double> in (3,5);
	in.Random ();

	std::cout << in << std::endl << std::endl;

	MXDump(in, "test.mat", "imat", "");

	Matrix<double> out;
	MXRead(out, "test.mat", "imat", "");

	std::cout << out << std::endl << std::endl;

	Matrix<cxfl> r1 (4,8);
	r1.Random ();
	MXDump (r1, "rtest.mat", "rmat", "");
	std::cout << r1 << std::endl << std::endl;
	
	Matrix<cxfl> r2;
	MXRead (r2, "rtest.mat", "rmat", "");
	std::cout << r2 << std::endl << std::endl;

#else

	std::cout << "MATLAB root not set during configuration (--with-matlabroot).\n Test skipped." << std::endl;

#endif

	return true;

}

