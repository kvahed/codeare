/*
 * t_hdf5.cpp
 *
 *  Created on: May 19, 2013
 *      Author: kvahed
 */

#include "MLFile.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

using namespace codeare::matrix::io;

std::string mname = "A";
std::string fname = "test.mat";

template <class T>
inline static bool write (const Matrix<T> A) {
	MLFile mfw (fname, WRITE);
	mfw.Write (A, mname);
}

template <class T>
inline static bool read (Matrix<T>& A) {
	MLFile mfr (fname, READ);
	A = mfr.Read<T>(mname);
}

template<class T>
inline static bool check () {

	Matrix<T> A = rand<T,uniform>(3,4), B;

	write(A);
	read(B);

	Matrix<unsigned short> C (A == B);

#if defined (VERBOSE)
	std::cout << C;
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


