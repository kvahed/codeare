/*
 * Linear.hpp
 *
 *  Created on: Apr 5, 2015
 *      Author: kvahed
 */

#ifndef _LINEAR_HPP_
#define _LINEAR_HPP_

#include "Matrix.hpp"
#include "Optimisable.hpp"

namespace codeare {
namespace optimisation {

template<class T>
class Linear {
public:
	Linear (Optimisable* opt = 0, const int& verbosity = 0) :
		_opt(opt), _verbosity(verbosity) {}
	virtual ~Linear () {}
	virtual int Minimise () = 0;
protected:
	int _verbosity;
	Optimisable* _opt;
};

}
}



#endif /* _LINEAR_HPP_ */
