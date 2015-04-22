/*
 * NLCG.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: kvahed
 */

#ifndef _NLCG_HPP_
#define _NLCG_HPP_

#include <NonLinear.hpp>

namespace codeare {
    namespace optimisation {
        
template<class T> class NLCG : public NonLinear<T> {

public:
	NLCG (const size_t& iterations) : NonLinear<T>::NonLinear (iterations)  {};
	virtual ~NLCG() {};
	inline Matrix<T> Minimize() { return Matrix<T>(); }
private:

};
        
    }}

#endif // _NLCG_HPP_
