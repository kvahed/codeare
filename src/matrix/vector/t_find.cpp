/*
 * t_find.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: kvahed
 */

#include "Creators.hpp"

template<class T> inline int check () {
	Vector<T> a = randn<T>(5,1).Container();
	std::cout << "a = [";
	std::cout << a << "]" << std::endl;
	std::cout << "find (a>0)" << std::endl;
	std::cout << find(a>0) << std::endl;
	std::cout << "find (a >-0.5 & a < 0.5)" << std::endl;
	std::cout << find(a >-0.5 & a < 0.5) << std::endl;
	return 0;
}

int main (int args, const char** argv) {
	return check<float>() + check<double>() + check<cxfl>() + check<cxdb>();
}
