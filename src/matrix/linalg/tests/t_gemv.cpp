#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

template<class T,class S> void gemv_check () {

    Matrix<T> A = rand<T> (4,4);
    Matrix<S> b = rand<S> (4,1);
    
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

    gemv_check<cxfl,cxfl>();
    gemv_check<cxdb,cxdb>();
    gemv_check<float,float>();
    gemv_check<double,double>();
    
    return 0;
    
}
