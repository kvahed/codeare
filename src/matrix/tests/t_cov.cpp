/*
 * t_cov.cpp
 *
 *  Created on: Oct 21, 2015
 *      Author: kvahed
 */

#include <Statistics.hpp>
#include <Print.hpp>
#include <Creators.hpp>
#include <Algos.hpp>

template<class T> inline int check () {
	Matrix<T> A = randn<T>(3,4);
    std::cout << "A = [" <<std::endl;
    std::cout << A << "];" << std::endl;
    std::cout << "cov(A)" <<std::endl;
    std::cout << cov(A) << std::endl;
	return 0;
}

int main (int narg, const char** argv) {
	return check<float>() + check<double>() + check<cxfl>() + check<cxdb>();
}


