/*
 * t_codeare.cpp
 *
 *  Created on: Jun 01, 2013
 *      Author: kvahed
 */

#include "DicomFile.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

#define VERBOSE
using namespace codeare::matrix::io;

std::string mname = "A";
std::string fname = "test.dcm";

template <class T> inline static bool write (const Matrix<T> A) {
	DicomFile dcm (fname, WRITE);
	dcm.Write (A, mname);
    return true;
}

template <class T> inline static bool read (Matrix<T>& A) {
	DicomFile mfr (fname, READ);
	A = mfr.Read<T>(mname);
    return true;
}

template<class T> inline static bool check () {

	Matrix<T> A = rand<T>(3,4), B;
	Matrix<cbool> C;

	write(A);
	read(B);

#if defined (VERBOSE)
	std::cout << A << std::endl;
	std::cout << B << std::endl;
#endif
	try {
		C = (A==B);
	} catch (const MatrixException&) {
		return false;
	}
#if defined (VERBOSE)
	std::cout << C << std::endl;
	std::cout << std::endl;
#endif

	return true;

}

int main (int args, char** argv) {

	if (check<float>() && check<double>() && check<cxfl>() && check<cxdb>())
        return 0;

	return 1;

}


