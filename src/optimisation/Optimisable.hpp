/*
 * Optimisable.hpp
 *
 *  Created on: Apr 5, 2015
 *      Author: kvahed
 */

#ifndef SRC_OPTIMISATION_OPTIMISABLE_HPP_
#define SRC_OPTIMISATION_OPTIMISABLE_HPP_

namespace codeare {
namespace optimisation {
class Optimisable {
public:
	virtual ~Optimisable () {};
	virtual int Objective () = 0;
};
}
}


#endif /* SRC_OPTIMISATION_OPTIMISABLE_HPP_ */
