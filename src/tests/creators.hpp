#include "matrix/Creators.hpp"

template <class T> bool
creatorstest (Connector<T>* rc) {

	Matrix<double> rd = rand<double> (8,3);
	std::cout << rd << std::endl;

	Matrix<float> rf = rand<float> (4,5);
	std::cout << rf << std::endl;

	Matrix<cxdb> cd = rand<cxdb> (3,7);
	std::cout << cd << std::endl;

	Matrix<cxfl> cf = rand<cxfl> (2,3);
	std::cout << cf << std::endl;

	Matrix<short> sh = rand<short> (5, 5);
	std::cout << sh << std::endl;

	Matrix<long> lg = rand<long> (6, 2);
	std::cout << lg << std::endl;

	return true;

}
