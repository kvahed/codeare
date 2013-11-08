#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

template<class T> void gemv_check () {

    Matrix<T> A = rand<T> (4,4);
    Matrix<T> b = rand<T> (4,1);
    
#ifdef VERBOSE
    std::cout << "A=\n" << A;
    std::cout << "b=\n" << b;
#endif
    
    A = gemv(A,b);

#ifdef VERBOSE
    std::cout << "A*b=\n" << A << std::endl;
#endif


}

int main (int args, char** argv) {

    gemv_check<cxfl>();
    gemv_check<cxdb>();
    gemv_check<float>();
    gemv_check<double>();
    
    return 0;
    
}
