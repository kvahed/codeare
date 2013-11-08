#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

template<class T> void sin_check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs = rand<T>(3,4);
    Matrix<T> lhs;

    lhs = sin (rhs);
    #ifdef VERBOSE
    std::cout << "A=\n" << rhs;
    std::cout << "sin(A)=\n" << lhs << std::endl;
    #endif    
}


int main (int args, char** argv) {

    sin_check<float>();
    sin_check<double>();
    sin_check<cxfl>();
    sin_check<cxdb>();
    
    return 0;
    
}
