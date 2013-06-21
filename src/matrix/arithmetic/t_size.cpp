#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Print.hpp"

template<class T> void size_check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs (3,4,5);
    Matrix<size_t> lhs;

    lhs = size(rhs);

#ifdef VERBOSE
    std::cout << "size(A)=\n" << lhs << std::endl;
#endif

}


int main (int args, char** argv) {

    size_check<float>();
    size_check<double>();
    
    return 0;
    
}
