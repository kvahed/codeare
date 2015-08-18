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
#include "Demangle.hpp"

namespace codeare {
    namespace optimisation {
        
    	template <class T>
        class NonLinear {
            
        public:
            NonLinear () {}
            NonLinear (const Params& params) {}
            NonLinear (const NonLinear& tocopy) {}
            virtual ~NonLinear () {}
            inline virtual void Minimise (Operator<T>* A, Matrix<T>& x) {}
            virtual std::ostream& Print (std::ostream& os) const {
                os << "  " << demangle(typeid(*this).name()).c_str() <<  std::endl;
                return os;
            };
            friend std::ostream& operator<< (std::ostream& os, const NonLinear<T>& oper) {
                return oper.Print(os);
            }
        protected:
            size_t _iterations;
            
        };
        
    }} // namespaces

#endif //_NON_LINEAR_
