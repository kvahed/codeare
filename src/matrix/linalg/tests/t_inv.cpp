#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

#define VERBOSE

template<class T> void inv_check () {

    Matrix<T> A = rand<T>(3,3);

#ifdef VERBOSE
    std::cout << "A=\n" << A << std::endl;
#endif
    
    A  = inv (A);
    
#ifdef VERBOSE
    std::cout << "inv(A)=\n" << A << std::endl;
    std::cout << std::endl;
#endif

}

int main (int args, char** argv) {

    inv_check<float>();
    inv_check<double>();
    inv_check<cxfl>();
    inv_check<cxdb>();
    
    return 0;
    
}
