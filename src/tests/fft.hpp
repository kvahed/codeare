#include "matrix/ft/DFT.hpp"

template <class T> bool
fftwtest (Connector<T>* rc) {

	typedef double Type;

	Matrix< std::complex<Type> > m   = phantom< std::complex<Type> > (128);
	Matrix< std::complex<Type> > k, i, j;
	Matrix<Type> mask = ones<Type> (size(m));
	Matrix<std::complex<Type> > pc = ones<std::complex<Type> > (size(m));

	DFT< Type > dft (size(m), mask, pc);

	k = dft   * m;
	i = dft ->* k;

	MXDump (m, "m.mat", "m");
	MXDump (k, "k.mat", "k");
	MXDump (i, "i.mat", "i");

	return true;

}

