#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

template<class T>
void emul_check () {

    Matrix<T> rhs = rand<T>(3,4);
    Matrix<T> lhs;

    lhs = rhs * rhs;
    #ifdef VERBOSE
    std::cout << "A=\n" << rhs;
    std::cout << "A.*A=\n" << lhs << std::endl;
    #endif
    
}


int main (int args, char** argv) {
    
    emul_check<float>();
    emul_check<double>();
    emul_check<cxfl>();
    emul_check<cxdb>();
    
    return 0;
    
}
