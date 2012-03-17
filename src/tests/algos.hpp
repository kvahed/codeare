#include "Statistics.hpp"
#include "Interpolate.hpp"

template <class T> bool
algotest (Connector<T>* rc) {

	Matrix<double> A;
	Matrix<cxdb>   B;
	Matrix<cxfl>   C;
	
	MXRead (A, std::string (base + std::string("real.mat")), "A");
	MXRead (B, std::string (base + std::string("complex.mat")), "B");

	Matrix<size_t> sz;

	sz = size(A);

	std::cout << "  size(A): \n    " << sz ;
	std::cout << "  numel(B): " << numel(B) << std::endl;
	
	C = (Matrix<cxfl>) Sum (B,4);
	sz = size(C);
	std::cout << "  size(C): \n    " << sz ;
	std::cout << "  numel(C): " << numel(C) << std::endl;

	MXDump (C, std::string (base + std::string("sumb3.mat")), "sumb3");

	C = (Matrix<cxfl>) SOS(B,4);
	MXDump (C, std::string (base + std::string("sumb3.mat")), "sos3");

	MXRead (B, std::string (base + std::string("covin.mat")), "A");
	B = inv(B = chol(B = cov(B)));
	MXDump (B, std::string (base + std::string("covout.mat")), "covA");

	Matrix<double> x  = Matrix<double>::LinSpace (0.0, 10.0, 8);
	Matrix<double> xi = Matrix<double>::LinSpace (0.0, 10.0, 100);

	Matrix<cxdb> Bi = interp1 (x, B, xi);

	MXDump (Bi, std::string (base + std::string("bi.mat")), "bi");

	return true;
	
}
