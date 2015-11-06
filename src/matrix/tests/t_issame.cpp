#include <Matrix.hpp>
#include <Creators.hpp>
#include <Symmetry.hpp>
#include <Lapack.hpp>
#include <Print.hpp>

template<class T> inline static int check () {
    Matrix<T> A = randn<T>(3,3);
    Matrix<T> B = A + 1e-12*randn<T>(3,3);
    Matrix<T> C = A(CR());
    int r1 = issame(A,C), r2 = issame(A,B);

    std::cout << "A = [" << std::endl;
    std::cout << A << "];" << std::endl;
    std::cout << "A(:):" << std::endl;
    std::cout << C << std::endl;
    std::cout << "issame(A,A(:)): " ;
    std::cout << r1 << std::endl << std::endl;
    std::cout << "issame(A,A+1e-12*randn(size(A))): " ;
    std::cout << r2 << std::endl << std::endl;

    if (typeid(T)==typeid(float) || typeid(T)==typeid(cxfl))
        return (r1==1&&r2==2);
    else
        return (r1==1&&r2==0);
}

int main (int args, char** argv) {
    if (check<float>() && check<double>() && check<cxfl>() && check<cxdb>())
        return 0;
    return 1;
}
