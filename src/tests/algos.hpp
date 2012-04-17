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
	
	C = (Matrix<cxfl>) sum (B,4);
	sz = size(C);
	std::cout << "  size(C): \n    " << sz ;
	std::cout << "  numel(C): " << numel(C) << std::endl;

	MXDump (C, std::string (base + std::string("sumb3.mat")), "sumb3");

	C = (Matrix<cxfl>) SOS(B,4);
	MXDump (C, std::string (base + std::string("sumb3.mat")), "sos3");

	MXRead (B, std::string (base + std::string("covin.mat")), "A");
	B = inv(B = chol(B = cov(B)));
	MXDump (B, std::string (base + std::string("covout.mat")), "covA");

	Matrix<double> x  = linspace<double> (0.0, 10.0, 8);
	Matrix<double> xi = linspace<double> (0.0, 10.0, 100);

	Matrix<cxdb> Bi = interp1 (x, B, xi, INTERP::AKIMA);

	MXDump (x, std::string (base + std::string("x.mat")), "x");
	MXDump (Bi, std::string (base + std::string("bi.mat")), "bi");
	MXDump (xi, std::string (base + std::string("xi.mat")), "xi");

	Matrix<double> y = linspace<double> (-0.2, 0.2, 5);
	Matrix<double> z = linspace<double> (-1.0,0.5,2);

	Matrix<double> xy = meshgrid<double> (x, y);
	Matrix<double> xyz = meshgrid<double> (x, y, z);

	MXDump (y, std::string (base + std::string("y.mat")), "y");
	MXDump (z, std::string (base + std::string("z.mat")), "z");
	MXDump (xy, std::string (base + std::string("xy.mat")), "xy");
	MXDump (xyz, std::string (base + std::string("xyz.mat")), "xyz");

	return true;
	
}
