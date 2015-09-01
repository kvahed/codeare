/*
 * SPLIT_BREGMAN.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: kvahed
 */

#ifndef _SPLIT_BREGMAN_HPP_
#define _SPLIT_BREGMAN_HPP_

#include <NonLinear.hpp>

template<class T> inline static Matrix<T>
shrink(const Matrix<T>& rhs, const typename TypeTraits<T>::RT& t) {
    Matrix<T> ret = rhs;
    typename TypeTraits<T>::RT rhsa = 0.;
    for (size_t i = 0; i < rhs.Size(); ++i) {
        rhsa = std::abs(rhs[i]);
        ret[i] = (rhsa > t) ? (1.0 - t/rhsa) * rhs[i] : T(0.); 
    }
    return ret;
}


namespace codeare {
    namespace optimisation {
        
        template<class T>
        class SplitBregman : public NonLinear<T> {
            
        public:
            SplitBregman () {};
            SplitBregman (const Params& p) : NonLinear<T>::NonLinear (p) {};
            SplitBregman (const SplitBregman& tocopy) {};            
            inline virtual void Minimise (Operator<T>* A, Matrix<T>& x) {};
            virtual ~SplitBregman() {};
        private:
            
        };
        
    }}

#endif // _SPLIT_BREGMAN_HPP_
