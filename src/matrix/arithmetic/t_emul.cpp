#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

template<class T>
void emul_check () {

    Matrix<T> A = rand<T>(3,4);
    Matrix<T> B = resize(transpose(A),size(A));

    B /= 2.0;
    B *= 2.0;
    B  = A * B;

#ifndef VERBOSE
    std::cout << "A=\n" << A;
    std::cout << "A.*A=\n" << B;
    std::cout << std::endl;
#endif
    
}


int main (int args, char** argv) {
    
    emul_check<float>();
    emul_check<double>();
    emul_check<cxfl>();
    emul_check<cxdb>();
    
    return 0;
    
}
