#include "Matrix.hpp"
#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "CX.hpp"
#include "Print.hpp"

template<class T, class S> void complex_check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs1 = rand<T>(3,4);
    Matrix<S> rhs2 = rand<S>(3,4);
    Matrix<std::complex<T> > lhs;

    lhs = complex(rhs1,rhs2);
    
    #ifdef VERBOSE
    std::cout << "A=\n" << rhs1;
    std::cout << "B=\n" << rhs2;
    std::cout << "A+i*B=\n" << lhs << std::endl;
    #endif
}


int main (int args, char** argv) {

    complex_check<float,float>();
    complex_check<float,double>();
    complex_check<double,float>();
    complex_check<double,double>();
    
    return 0;
    
}
