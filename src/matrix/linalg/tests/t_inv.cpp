#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

template<class T> void inv_check () {

    Matrix<T> A = rand<T,uniform>(3,3);

#ifdef VERBOSE
    std::cout << "A=\n" << A;
#endif
    
    A  = inv (A);
    
#ifdef VERBOSE
    std::cout << "inv(A)=\n" << A;
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
