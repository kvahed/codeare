/*
 * NonLinear.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: kvahed
 */

#ifndef _NON_LINEAR_
#define _NON_LINEAR_

#include "Matrix.hpp"
#include "Operator.hpp"
#include "Params.hpp"
#include "Lapack.hpp"

namespace codeare {
    namespace optimisation {
        
    	template <class T>
        class NonLinear {
            
        public:
            NonLinear (const Params& p) {}
            virtual ~NonLinear () {}
            virtual void Minimise (Operator<T>& A, Matrix<T>& x) {}
        protected:
            size_t _iterations;
            
        };
        
    }} // namespaces

#endif //_NON_LINEAR_
