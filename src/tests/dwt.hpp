#include "matrix/DWT.hpp"

template <class T> bool
dwttest (Connector<T>* rc) {

	Matrix<cxfl> m   = Matrix<cxfl>::Phantom2D(512);
	Matrix<cxfl> k, i;

	DWT dwt (m.Height(), WL_DAUBECHIES);

	k = dwt.Trafo(m);
	i = dwt.Adjoint(k);

#ifdef HAVE_MAT_H	

	MATFile* mf = matOpen ("dwtout.mat", "w");
	
	if (mf == NULL) {
		printf ("Error creating file %s\n", "");
		return false;
	}

	IO::MXDump (m, mf, "m", "");
	IO::MXDump (k, mf, "k", "");
	IO::MXDump (i, mf, "i", "");
	
	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", "");
		return false;
	}

#endif
	

	return true;

}

