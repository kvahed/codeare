/*
 * LBFGS.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: kvahed
 */

#ifndef _LBFGS_HPP_
#define _LBFGS_HPP_

#include <NonLinear.hpp>

namespace codeare {
    namespace optimisation {
        
        template<class T>
        class LBFGS : public NonLinear<T> {
            
        public:
            LBFGS (size_t iterations) : NonLinear<T>::NonLinear (iterations)  {};
            virtual ~LBFGS() {};
        private:
            
        };
        
    }}

#endif // _LBFGS_HPP_
