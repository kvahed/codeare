#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Statistics.hpp"

template<class T>
void emul_check () {

    Matrix<T> A = rand<T,uniform>(3,4);
    Matrix<T> B = rand<T,uniform>(3,4);

#ifdef VERBOSE
    std::cout << "A=\n" << A;
    std::cout << "B=\n" << B;
    std::cout << "nrmse(A,-B)=" << nrmse(A,-B);
    std::cout << std::endl;
#endif
    
}


int main (int args, char** argv) {
    
    emul_check<float>();
    emul_check<double>();
    emul_check<cxfl>();
    emul_check<cxdb>();
    
    return 0;
    
}
