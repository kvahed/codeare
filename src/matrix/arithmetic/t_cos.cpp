#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

template<class T>
void cos_check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs = rand<T,uniform>(3,4);
    Matrix<T> lhs;

    lhs = cos (rhs);
    #ifdef VERBOSE
    std::cout << "A=\n" << rhs;
    std::cout << "cos(A)=\n" << lhs << std::endl;
    #endif
    
}


int main (int args, char** argv) {
    
    cos_check<float>();
    cos_check<double>();
    cos_check<cxfl>();
    cos_check<cxdb>();
    
    return 0;
    
}
