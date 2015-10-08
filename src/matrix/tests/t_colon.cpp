#include "Matrix.hpp"
#include "Creators.hpp"
#include "Print.hpp"

typedef Range<false> R;
typedef Range<true> CR;

int main (int narg, char** argv) {
    Matrix<float> A = randn<float>(3,4);

    std::cout << __LINE__ << std::endl;
    std::cout << "A" << std::endl;
    std::cout << A << std::endl;
    std::cout << __LINE__ << std::endl;
    std::cout << "A(Range(0,2),Range(2,2))" << std::endl;
    std::cout << Matrix<float>(A(CR(0,2),CR(2,2))) << std::endl;
    std::cout << __LINE__ << std::endl;
    std::cout << "A(Range(0,2),Range(2))" << std::endl;
    std::cout << Matrix<float>(A(CR(0,2),CR(2))) << std::endl;
    std::cout << __LINE__ << std::endl;
    return 0;
}
