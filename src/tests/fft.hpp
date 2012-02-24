#include "matrix/FFT.hpp"
#include "matrix/DFT.hpp"

template <class T> bool
fftwtest (Connector<T>* rc) {

	std::string in   = std::string (base + std::string ("/in.mat"));
	std::string iout = std::string (base + std::string ("/iout.mat"));
	std::string kout = std::string (base + std::string ("/kout.mat"));
	
	Matrix<cxfl> m   = Matrix<cxfl>::Phantom2D(512);
	Matrix<cxfl> k, i, j;

	k = FFT::Forward(m);
	i = FFT::Backward(k);

	Matrix<double> msk;
	Matrix<double> pdf;
	Matrix<cxfl>   dat;
	Matrix<cxfl>   phc;
	Matrix<cxfl>   tst;

	IO::MXRead (msk, "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "mask");
	IO::MXRead (pdf,  "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "pdf");
	IO::MXRead (dat, "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "data");
	IO::MXRead (phc, "/Users/kvahed/git/codeare/share/compressedsensing/noisydata.mat", "ph");


	dat /= pdf;
	DFT dft (2, 256, msk, phc);

	dat = dft.Adjoint (dat);
	tst = dft.Trafo   (dat);
	tst = dft.Adjoint (tst);
	tst = dft.Trafo   (tst);

#ifdef HAVE_MAT_H	

	MATFile* mf = matOpen ("fftout.mat", "w");
	
	if (mf == NULL) {
		printf ("Error creating file %s\n", "");
		return false;
	}

	IO::MXDump (m, mf, "m", "");
	IO::MXDump (k, mf, "k", "");
	IO::MXDump (i, mf, "i", "");
	IO::MXDump (dat, mf, "dat", "");
	IO::MXDump (tst, mf, "tst", "");
	
	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", "");
		return false;
	}

#endif
	

	return true;

}

