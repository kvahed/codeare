#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

template<class T> void tan_check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs = rand<T>(3,4);
    Matrix<T> lhs;

    lhs = tan (rhs);
    #ifdef VERBOSE
    std::cout << "A=\n" << rhs;
    std::cout << "tan(A)=\n" << lhs << std::endl;
    #endif
}


int main (int args, char** argv) {

    tan_check<float>();
    tan_check<double>();
    tan_check<cxfl>();
    tan_check<cxdb>();
    
    return 0;
    
}
