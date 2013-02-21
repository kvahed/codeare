#include "Matrix.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

template<class T,class S> void gemm_check () {

    Matrix<T> A = rand<T> (4,3);
    Matrix<S> B = rand<S> (3,4);
    
#ifdef VERBOSE
    std::cout << "A=\n" << A;
    std::cout << "B=\n" << B;
#endif
    
    A = gemm(A,B);

#ifdef VERBOSE
    std::cout << "A*B=\n" << A << std::endl;
#endif

    A = gemm(A,B,'T','T');

#ifdef VERBOSE
    std::cout << "(A*B)'*B'=\n" << A << std::endl;
#endif

}

int main (int args, char** argv) {

    gemm_check<cxfl,cxfl>();
    gemm_check<cxdb,cxdb>();
    gemm_check<float,float>();
    gemm_check<double,double>();
    
    return 0;
    
}
