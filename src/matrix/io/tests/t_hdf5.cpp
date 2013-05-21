/*
 * t_hdf5.cpp
 *
 *  Created on: May 19, 2013
 *      Author: kvahed
 */

#include "HDF5File.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

using namespace codeare::matrix::io;

std::string mname = "/group1/group2/A";
std::string fname = "test.h5";

template <class T>
inline static bool write (const Matrix<T> A) {
	HDF5File h5fw (fname, WRITE);
	h5fw.Write (A, mname);
}

template <class T>
inline static bool read (Matrix<T>& A) {
	HDF5File h5fr (fname, READ);
	A = h5fr.Read<T>(mname);
}

template<class T>
inline static bool check () {

	Matrix<T> A = rand<T>(3,4), B;

	write(A);
	read(B);

#if defined (VERBOSE)
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

	return 0;

}


