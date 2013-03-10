#include "Matrix.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

template<class T> void dot_check () {

    Matrix<T> A = rand<T> (4,1);
    Matrix<T> B = rand<T> (4,1);
    
#ifdef VERBOSE
    std::cout << "A=\n" << A;
    std::cout << "B=\n" << B;
    std::cout << "A * B=\n" << dot(A,B) << std::endl;
    std::cout << "A**H * B=\n" << dotc(A,B) << std::endl;
#endif

}

int main (int args, char** argv) {

    dot_check<cxfl>();
    dot_check<cxdb>();
    dot_check<float>();
    dot_check<double>();
    
    return 0;
    
}
