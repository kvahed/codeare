#include "Matrix.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

template<class T> void pinv_check () {

    Matrix<T> A = rand<T>(4,3);

#ifdef VERBOSE
    std::cout << "A=\n" << A;
#endif
    
    A = pinv (A);
    
#ifdef VERBOSE
    std::cout << "pinv(A)=\n" << A << std::endl;
#endif
    
}

int main (int args, char** argv) {

    pinv_check<cxfl>();
    pinv_check<cxdb>();
    pinv_check<float>();
    pinv_check<double>();
    
    return 0;
    
}
