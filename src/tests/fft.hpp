#include "matrix/FFT.hpp"
#include "matrix/DFT.hpp"

template <class T> bool
fftwtest (Connector<T>* rc) {

	std::string in   = std::string (base + std::string ("/in.mat"));
	std::string iout = std::string (base + std::string ("/iout.mat"));
	std::string kout = std::string (base + std::string ("/kout.mat"));
	
	Matrix<cxfl> m   = Matrix<cxfl>::Phantom2D(512);
	Matrix<cxfl> k, i, j;

	ticks cgstart = getticks();
	for (size_t n = 0; n < 50; n++)
		k = FFT::Forward(m);
	printf ("FFT: %.4f s\n", elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());
	i = FFT::Backward(k);

	Matrix<int> size (2,1);
	size[0] = 512;
	size[1] = 512;

	Matrix<double> mask;

	IO::MXRead (mask, "/Users/kvahed/git/codeare/share/compressedsensing/brain512.mat", "mask");

	DFT dft (size, mask);

	cgstart = getticks();
	for (size_t n = 0; n < 50; n++)
		j = dft.Trafo(m);
	printf ("DFT: %.4f s\n", elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());


#ifdef HAVE_MAT_H	

	MATFile* mf = matOpen ("fftout.mat", "w");
	
	if (mf == NULL) {
		printf ("Error creating file %s\n", "");
		return false;
	}

	IO::MXDump (m, mf, "m", "");
	IO::MXDump (k, mf, "k", "");
	IO::MXDump (i, mf, "i", "");
	IO::MXDump (j, mf, "j", "");
	
	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", "");
		return false;
	}

#endif
	

	return true;

}

