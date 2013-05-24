/*
 * t_hdf5.cpp
 *
 *  Created on: May 19, 2013
 *      Author: kvahed
 */

#include "NIFile.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

using namespace codeare::matrix::io;

std::string fname = "test.nii";
std::string mname = "t_nifti";

template <class T>
inline static bool write (const Matrix<T> A) {
	NIFile nf (fname, WRITE);
	nf.Write (A, mname);
}

template <class T>
inline static bool read (Matrix<T>& A) {
	NIFile nf (fname, READ);
	A = nf.Read<T>(mname);
}

template<class T>
inline static bool check () {

	Matrix<T> A = rand<T>(3,4), B;

	write(A);
	read(B);

/*	h5write (A,fname);
    A = h5read<T> (fname);*/

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


