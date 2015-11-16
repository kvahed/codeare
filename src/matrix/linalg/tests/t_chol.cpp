#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"
#include "Print.hpp"

#define VERBOSE

template<class T> void chol_check () {

    Matrix<T> A = rand<T> (4,3);
    
#ifdef VERBOSE
    std::cout << "A=\n" << A << std::endl;
#endif
    
    A = A.mult(A); // m*m' (Must be positive definite)
    A = chol (A);

#ifdef VERBOSE
    std::cout << "chol(A)=\n" << A  << std::endl;
#endif

}

int main (int args, char** argv) {

    chol_check<cxfl>();
    chol_check<cxdb>();
    chol_check<float>();
    chol_check<double>();
    
    return 0;
    
}
