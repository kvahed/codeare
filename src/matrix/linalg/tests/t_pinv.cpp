#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"
#include "Print.hpp"

template<class T> void check () {

    Matrix<T> A = rand<T>(8,3);
    
#ifdef VERBOSE
    std::cout << "A=\n" << A << std::endl;
#endif
    
    Matrix<T> B = pinv (A);
    
#ifdef VERBOSE
    std::cout << "pinv(A)=\n" << B << std::endl;
    std::cout << std::endl;
#endif
    
}

int main (int args, char** argv) {

    check<float>();
    check<double>();
    check<cxfl>();
    check<cxdb>();

    
    return 0;
    
}
