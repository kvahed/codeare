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
    return true;
}

template <class T>
inline static bool read (Matrix<T>& A) {
	MLFile mfr (fname, READ);
	A = mfr.Read<T>(mname);
    return true;
}

template<class T>
inline static bool check () {

	Matrix<T> A = rand<T>(3,4), B;
	write(A);
	read(B);
	return issame(A,B);

}

int main (int args, char** argv) {
    if (args == 1) {
        if (check<float>() && check<double>() && check<cxfl>() && check<cxdb>())
            return 0;
    } else {
        MLFile mf (argv[1],READ);
        mf.Read();
        std::cout << wspace << std::endl;
        return 0;
    }
    return 1;
}


