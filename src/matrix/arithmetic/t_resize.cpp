#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Print.hpp"

template<class T>
void check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs = rand<T,uniform>(3,4);
    Matrix<T> blhs, slhs;

    slhs = resize (rhs,3,3);
    blhs = resize (rhs,7,2);
    #ifdef VERBOSE
    std::cout << "A=\n" << rhs;
    std::cout << "resize(A,3,3)=\n" << slhs;
    std::cout << "resize(A,2,5)=\n" << blhs << std::endl;
    #endif
    
}


int main (int args, char** argv) {
    
    check<float>();
    check<double>();
    check<cxfl>();
    check<cxdb>();
    check<short>();
    check<long>();

    return 0;
    
}
