/*
 * t_codeare.cpp
 *
 *  Created on: Jun 01, 2013
 *      Author: kvahed
 */

#include "CODFile.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

using namespace codeare::matrix::io;

std::string mname = "A";
std::string fname = "test.cod";

template <class T>
inline static bool write (const Matrix<T> A) {
	CODFile mfw (fname, WRITE);
	mfw.Write (A, mname);
}

template <class T>
inline static bool read (Matrix<T>& A) {
	CODFile mfr (fname, READ);
	A = mfr.Read<T>(mname);
}

template<class T>
inline static bool check () {

	Matrix<T> A = rand<T>(3,4), B;

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

	if (check<float>() && check<double>() && check<cxfl>() && check<cxdb>())
        return 0;

	return 1;

}


