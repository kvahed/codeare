#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "CX.hpp"

template<class T> void abs_check () {

    Matrix<T> rhs = rand<T>(3,4);
    Matrix<T> lhs;

    lhs = abs (rhs);
    #ifdef VERBOSE
    std::cout << "A=\n" << rhs;
    std::cout << "abs(A)=\n" << lhs << std::endl;
    #endif
}


int main (int args, char** argv) {

    abs_check<float>();
    abs_check<double>();
    abs_check<cxfl>();
    abs_check<cxdb>();
    
    return 0;
    
}
