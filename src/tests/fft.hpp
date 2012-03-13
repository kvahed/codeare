#include "matrix/FFT.hpp"

template <class T> bool
fftwtest (Connector<T>* rc) {

	Matrix<cxfl> m   = Matrix<cxfl>::Phantom2D(512);
	Matrix<cxfl> k, i, j;
	DFT<cxfl> dft (size(m));
	
	for (size_t l = 0; l < 100; l++)
		k = dft * m;

	for (size_t l = 0; l < 100; l++)
		i = dft ->* k;

	/*
	Matrix<float> msk;
	Matrix<float> pdf;
	Matrix<cxfl>   dat;
	Matrix<cxfl>   phc;
	Matrix<cxfl>   tst;

	MXRead (msk, "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "mask");
	MXRead (pdf,  "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "pdf");
	MXRead (dat, "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "data");
	MXRead (phc, "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "ph");

	dat /= pdf;
	DFT<cxfl> dft (2, 256, msk, phc);

	dat = dft * dat;
	tst = dft * dat;
	tst = dft->*tst;
	tst = dft * tst;
	*/

#ifdef HAVE_MAT_H	

	MATFile* mf = matOpen ("fftout.mat", "w");
	
	if (mf == NULL) {
		printf ("Error creating file %s\n", "");
		return false;
	}

	MXDump (m, mf, "m", "");
	MXDump (k, mf, "k", "");
	MXDump (i, mf, "i", "");
	/*
	MXDump (dat, mf, "dat", "");
	MXDump (tst, mf, "tst", "");
	*/
	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", "");
		return false;
	}

#endif
	

	return true;

}

