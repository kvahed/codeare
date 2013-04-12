#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

template<class T> void inv_check () {

    Matrix<T> A = rand<T>(4,4);

#ifdef VERBOSE
    std::cout << "A=\n" << A;
#endif
    
    A  = inv (A);
    
#ifdef VERBOSE
    std::cout << "inv(A)=\n" << A << std::endl;
#endif

}

int main (int args, char** argv) {

    inv_check<cxfl>();
    inv_check<cxdb>();
    inv_check<float>();
    inv_check<double>();
    
    return 0;
    
}
