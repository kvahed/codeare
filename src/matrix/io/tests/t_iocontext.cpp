/*
 * t_hdf5.cpp
 *
 *  Created on: May 19, 2013
 *      Author: kvahed
 */

#include "Algos.hpp"
#include "Creators.hpp"
#include "IOContext.hpp"

using namespace codeare::matrix::io;

std::string mname = "/group1/group2/A";
std::string fname = "test.h5";

template<class T>
inline static bool check () {

	Matrix<T> A = rand<T,uniform>(3,4), B;

    IOContext ioc (fname);

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


