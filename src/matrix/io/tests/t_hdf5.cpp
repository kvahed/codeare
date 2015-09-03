/*
 * t_hdf5.cpp
 *
 *  Created on: May 19, 2013
 *      Author: kvahed
 */

#include "Matrix.hpp"
#include "HDF5File.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Print.hpp"

using namespace codeare::matrix::io;

std::string mname = "A";
std::string fname = "test.h5";

template <class T>
inline static bool write (const Matrix<T> A) {
	HDF5File nf (fname, WRITE);
	nf.Write (A, mname);
    return true;
}

template <class T>
inline static bool read (Matrix<T>& A) {
	HDF5File nf (fname, READ);
	A = nf.Read<T>(mname);
    return true;
}

template<class T>
inline static bool check () {

	Matrix<T> A = rand<T>(3,4), B;

	// Class interface
	write(A);
	read(B);

#if defined (VERBOSE)
	std::cout << A;
	std::cout << B;
	std::cout << (A == B);
	std::cout << std::endl;
#endif

	// Convenience interface
	h5write (A,fname);
    B = h5read<T> (fname,"/A");

#if defined (VERBOSE)
	std::cout << A;
	std::cout << B;
	std::cout << (A == B);
	std::cout << std::endl;
#endif

	return true;

}

int main (int args, char** argv) {

	if (!check<float>())
		return 1;
	if (!check<double>())
		return 1;
	if (!check<cxfl>())
		return 1;
	if (!check<cxdb>())
		return 1;
	if (!check<short>())
		return 1;

	return 0;

}


