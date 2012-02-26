#include "matrix/FFT.hpp"
#include "matrix/DFT.hpp"

template <class T> bool
fftwtest (Connector<T>* rc) {

	std::string in   = std::string (base + std::string ("/in.mat"));
	std::string iout = std::string (base + std::string ("/iout.mat"));
	std::string kout = std::string (base + std::string ("/kout.mat"));
	
	Matrix<cxfl> m   = Matrix<cxfl>::Phantom2D(512);
	Matrix<cxfl> k, i, j;

	k = fft (m);
	i = ifft (k);

	Matrix<double> msk;
	Matrix<double> pdf;
	Matrix<cxfl>   dat;
	Matrix<cxfl>   phc;
	Matrix<cxfl>   tst;

	MXRead (msk, "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "mask");
	MXRead (pdf,  "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "pdf");
	MXRead (dat, "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "data");
	MXRead (phc, "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "ph");


	dat /= pdf;
	DFT dft (2, 256, msk, phc);

	//dat = dft.Trafo (dat);
	dat = dft * dat;
	tst = dft.Trafo (dat);
	tst = dft.Adjoint (tst);
	tst = dft.Trafo (tst);

#ifdef HAVE_MAT_H	

	MATFile* mf = matOpen ("fftout.mat", "w");
	
	if (mf == NULL) {
		printf ("Error creating file %s\n", "");
		return false;
	}

	MXDump (m, mf, "m", "");
	MXDump (k, mf, "k", "");
	MXDump (i, mf, "i", "");
	MXDump (dat, mf, "dat", "");
	MXDump (tst, mf, "tst", "");
	
	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", "");
		return false;
	}

#endif
	

	return true;

}

