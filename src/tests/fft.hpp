#include "matrix/ft/DFT.hpp"

template <class T> bool
fftwtest (Connector<T>* rc) {

	Matrix< std::complex<float> > m   = phantom< std::complex<float> > (256);
	Matrix< std::complex<float> > k, i, j;
	DFT< float > dft (size(m));

	k = dft   * m;
	i = dft ->* k;

	MXDump (m, "m.mat", "m");
	MXDump (k, "k.mat", "k");
	MXDump (i, "i.mat", "i");

	return true;

}

