/*
 * SPLIT_BREGMAN.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: kvahed
 */

#ifndef _SPLIT_BREGMAN_HPP_
#define _SPLIT_BREGMAN_HPP_

#include <NonLinear.hpp>

namespace codeare {
    namespace optimisation {
        
        template<class T>
        class SplitBregman : public NonLinear<T> {
            
        public:
            SplitBregman (const size_t& iterations) : NonLinear<T>::NonLinear (iterations) {};
            virtual ~SplitBregman() {};
        private:
            
        };
        
    }}

#endif // _SPLIT_BREGMAN_HPP_
