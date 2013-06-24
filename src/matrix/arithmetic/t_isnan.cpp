#include "Trigonometry.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Print.hpp"

template<class T> void isnan_check () {

    using namespace codeare::matrix::arithmetic;
    
    Matrix<T> rhs (3,2);
    Matrix<cbool> lhs;

    rhs(1,0) = 1.0/0.0;
    rhs(0,1) = 1.0/0.0;
    lhs = isnan(rhs);

#ifdef VERBOSE
    std::cout << "size(A)=\n" << lhs << std::endl;
#endif

}


int main (int args, char** argv) {

    isnan_check<float>();
    isnan_check<double>();
    isnan_check<cxfl>();
    isnan_check<cxdb>();
    /*   isnan_check<short>();
    isnan_check<long>();
    isnan_check<unsigned short>();
    isnan_check<size_t>();
    */
    return 0;
    
}
