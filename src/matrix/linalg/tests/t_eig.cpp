#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"
#include "Print.hpp"

#define VERBOSE

template<class T> void eig_check () {

    Matrix<T> A = rand<T>(4,4);
    eig_t<T> e;

#ifdef VERBOSE    
    std::cout << "A =\n" << A << std::endl;
#endif
    
    e = eig2<T> (A, 'V', 'V');
    
#ifdef VERBOSE    
    std::cout << "ev=\n" << e.ev << std::endl;
    std::cout << "ev=\n" << eig(A) << std::endl;
    std::cout << "lv=\n" << e.lv << std::endl;
    std::cout << "rv=\n" << e.rv << std::endl;
    std::cout << std::endl;
#endif

    A = gemm(A,A,'N','C');

#ifdef VERBOSE    
    std::cout << "A*A' =\n" << A << std::endl;
#endif
    
    e = eig2<T> (A, 'V', 'V');
    
#ifdef VERBOSE    
    std::cout << "ev=\n" << e.ev << std::endl;
    std::cout << "ev=\n" << eig(A) << std::endl;
    std::cout << "lv=\n" << e.lv << std::endl;
    std::cout << "rv=\n" << e.rv << std::endl;
    std::cout << std::endl;
#endif

}

int main (int args, char** argv) {

    eig_check<float>();
    eig_check<double>();
    eig_check<cxfl>();
    eig_check<cxdb>();
    
    return 0;
    
}
