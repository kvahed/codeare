template <class T> bool
iotest (Connector<T>* rc) {

	Matrix<double> d (100,300);
	d.Random();

	Matrix<cxfl> c (100,300);
	c.Random();

	PRDump (d, "d.cod");
	PRDump (c, "c.cod");

	PRRead (c, "c.cod");

#ifdef HAVE_MAT_H	
	MXDump (d, "test.mat", "d");
	MXDump (c, "test.mat", "c");
#endif

	return true;

}
