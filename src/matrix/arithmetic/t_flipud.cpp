#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

template<class T>
void check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs = rand<T,uniform>(3,4);
    Matrix<T> lhs;

    lhs = flipud (rhs);
    #ifdef VERBOSE
    std::cout << "A=\n" << rhs;
    std::cout << "flipud(A)=\n" << lhs << std::endl;
    #endif
    
}


int main (int args, char** argv) {
    
    check<float>();
    check<double>();
    check<cxfl>();
    check<cxdb>();
    check<short>();
    check<long>();
    
    return 0;
    
}
