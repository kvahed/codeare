#include "matrix/FFT.hpp"

template <class T> bool
fftwtest (Connector<T>* rc) {

	std::string in   = std::string (base + std::string ("/in.mat"));
	std::string iout = std::string (base + std::string ("/iout.mat"));
	std::string kout = std::string (base + std::string ("/kout.mat"));
	
	Matrix<cxfl> m   = Matrix<cxfl>::Phantom2D(512);
	Matrix<cxfl> k, i;

	k = FFT::Forward(m);
	i = FFT::Backward(k);

#ifdef HAVE_MAT_H	

	MATFile* mf = matOpen ("fftout.mat", "w");
	
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

