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

template <class T> inline static bool write (const Matrix<T> A) {
	HDF5File nf (fname, WRITE);
	nf.Write (A, mname);
    return true;
}

template <class T> inline static bool read (Matrix<T>& A) {
	HDF5File nf (fname, READ);
	A = nf.Read<T>(mname);
    return true;
}

template<class T> inline static bool check () {

	Matrix<T> A = rand<T>(3,4), B;

	// Class interface
	write(A);
	read(B);

	std::cout << (A == B);
	std::cout << std::endl;

	// Convenience interface
	h5write (A,fname);
    B = h5read<T> (fname,"/A");

 	std::cout << (A == B);
	std::cout << std::endl;

	return true;

}

int main (int args, char** argv) {

    if (args == 1) {
        if (check<float>() && check<double>() && check<cxfl>() && check<cxdb>())
            return 0;
    } else {
        HDF5File h5f (argv[1]);
        h5f.Read();
        std::cout << wspace << std::endl;
        return 0;
    }
    return 1;

}


