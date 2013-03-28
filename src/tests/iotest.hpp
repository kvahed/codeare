#include "matrix/Creators.hpp"

template <class T> bool
iotest (Connector<T>* rc) {

	Matrix<double> d = rand<double> (100,300);

	Matrix<cxfl> c = rand<cxfl> (100,300);
	c.SetClassName ("ComlexSingle");

	PRDump (d, "d.cod");
	PRDump (c, "c.cod");

	PRRead (c, "c.cod");

#ifdef HAVE_MAT_H	
	MXDump (d, "test.mat", "d");
	MXDump (c, "test.mat", "c");
#endif

	return true;

}
