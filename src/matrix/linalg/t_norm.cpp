#include "Matrix.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

template<class T> void norm_check () {

    Matrix<T> A = rand<T> (4,3);
    
#ifdef VERBOSE
    std::cout << "A=\n" << A;
#endif
    

#ifdef VERBOSE
    std::cout << "norm(A)=" << norm(A) << std::endl;
#endif


}

int main (int args, char** argv) {

    norm_check<cxfl>();
    norm_check<cxdb>();
    norm_check<float>();
    norm_check<double>();
    
    return 0;
    
}
