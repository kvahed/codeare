#include "matrix/DWT.hpp"

template <class T> bool
dwttest (Connector<T>* rc) {

	Matrix<cxfl> m   = Matrix<cxfl>::Phantom2D(512);
	Matrix<cxfl> k, i;

	k = DWT::Forward(m);
	i = DWT::Backward(k);

#ifdef HAVE_MAT_H	

	MATFile* mf = matOpen ("dwtout.mat", "w");
	
	if (mf == NULL) {
		printf ("Error creating file %s\n", "");
		return false;
	}

	m.MXDump (mf, "m", "");
	k.MXDump (mf, "k", "");
	i.MXDump (mf, "i", "");
	
	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", "");
		return false;
	}

#endif
	

	return true;

}

