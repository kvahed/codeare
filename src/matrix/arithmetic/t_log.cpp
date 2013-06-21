#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

template<class T> void log_check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs = rand<T,uniform>(3,4);
    Matrix<T> lhs;

    lhs = log (rhs);
    #ifdef VERBOSE
    std::cout << "A=\n" << rhs;
    std::cout << "log(A)=\n" << lhs << std::endl;
    #endif
}


int main (int args, char** argv) {

    log_check<float>();
    log_check<double>();
    log_check<cxfl>();
    log_check<cxdb>();
    
    return 0;
    
}
