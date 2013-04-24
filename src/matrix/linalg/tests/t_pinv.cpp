#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

template<class T> void pinv_check () {

    Matrix<T> A = rand<T>(8,3);
    
#ifdef VERBOSE
    std::cout << "A=\n" << A;
#endif
    
    Matrix<T> B = pinv (A);
    
#ifdef VERBOSE
    std::cout << "pinv(A)=\n" << B << std::endl;
#endif
    
}

int main (int args, char** argv) {

    pinv_check<cxfl>();
    pinv_check<cxdb>();
    pinv_check<float>();
    pinv_check<double>();
    
    return 0;
    
}
