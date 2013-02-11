#include "Matrix.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

template<class T, class S> void svd_check () {

    Matrix<T> A = rand<T>(3,4), U, V;
    Matrix<S> s;

    #ifdef VERBOSE
    std::cout << "A=\n" << A;
    #endif
    fflush (stdout);
    
    svd (A, s, U, V, 'A');
    
    #ifdef VERBOSE
    std::cout << "U=\n" << U;
    std::cout << "S=\n" << s;
    std::cout << "V=\n" << V << std::endl;
    #endif

}

int main (int args, char** argv) {

    svd_check<cxfl,float>();
    svd_check<cxdb,double>();
    svd_check<float,float>();
    svd_check<double,double>();
    
    return 0;
    
}
