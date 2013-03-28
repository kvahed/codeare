#include "Trigonometry.hpp"
#include "Creators.hpp"

template<class T> void log10_check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs = rand<T>(3,4);
    Matrix<T> lhs;

    lhs = log10 (rhs);
    #ifdef VERBOSE
    std::cout << "A=\n" << rhs;
    std::cout << "log10(A)=\n" << lhs << std::endl;
    #endif    
}


int main (int args, char** argv) {

    log10_check<float>();
    log10_check<double>();
    log10_check<cxfl>();
    log10_check<cxdb>();
    
    return 0;
    
}
