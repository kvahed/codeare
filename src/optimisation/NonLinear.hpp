/*
 * NonLinear.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: kvahed
 */

#ifndef _NON_LINEAR_
#define _NON_LINEAR_

#include "Matrix.hpp"

namespace codeare {
    namespace optimisation {
        
    	template <class T>
        class NonLinear {
            
        public:
            NonLinear (const size_t& iterations) : _iterations(iterations) {}
            virtual ~NonLinear () {}
        protected:
            size_t _iterations;
            
        };
        
    }} // namespaces

#endif //_NON_LINEAR_
