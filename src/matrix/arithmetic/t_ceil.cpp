#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

template<class T> void exp_check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs = 10. * ( rand<T,uniform>(3,4) + 1.0);
    Matrix<size_t> lhs;

    lhs = (Matrix<size_t>) ceil (rhs);

#ifdef VERBOSE
    std::cout << "A=\n" << rhs;
    std::cout << "ceil(A)=\n" << lhs << std::endl;
#endif

}


int main (int args, char** argv) {

    exp_check<float>();
    exp_check<double>();
    
    return 0;
    
}
