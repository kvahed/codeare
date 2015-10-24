#include <Matrix.hpp>
#include <Creators.hpp>
#include <Symmetry.hpp>
#include <Lapack.hpp>
#include <Print.hpp>

template<class T> inline static int check () {
    Matrix<T> A = randn<T>(3,3);
    std::cout << "A = [" << std::endl;
    std::cout << A << "];" << std::endl;
    std::cout << "issymmtric(A): " ;
    std::cout << issymmetric(A) << std::endl << std::endl;
    A = gemm(A,A,'N','C');
    std::cout << "A*A' = [" << std::endl;
    std::cout << A << "];" << std::endl;
    std::cout << "issymmtric(A*A'): " ;
    std::cout << issymmetric(A) << std::endl << std::endl;
    return 0;
}

int main (int args, char** argv) {
    return check<float>() + check<double>();
}
