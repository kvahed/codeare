/*
 * Linear.hpp
 *
 *  Created on: Apr 5, 2015
 *      Author: kvahed
 */

#ifndef _LINEAR_HPP_
#define _LINEAR_HPP_

#include "Matrix.hpp"
#include "Operator.hpp"

namespace codeare {
namespace optimisation {

template<class T>
class Linear {
public:
	Linear (const int& verbosity = 0) :	_verbosity(verbosity) {}
	virtual ~Linear () {}
	virtual Matrix<T> Solve(const Operator<T>& A, const MatrixType<T>& b) { return Matrix<T>(); }

protected:
	int _verbosity;
	//const Operator<T>& _E;
};

}
}



#endif /* _LINEAR_HPP_ */
