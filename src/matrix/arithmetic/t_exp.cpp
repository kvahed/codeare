#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

template<class T> void exp_check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs = rand<T,uniform>(3,4);
    Matrix<T> lhs;

    lhs = exp (rhs);
    #ifdef VERBOSE
    std::cout << "A=\n" << rhs;
    std::cout << "exp(A)=\n" << lhs << std::endl;
    #endif
}


int main (int args, char** argv) {

    exp_check<float>();
    exp_check<double>();
    exp_check<cxfl>();
    exp_check<cxdb>();
    
    return 0;
    
}
