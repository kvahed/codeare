#include <Matrix.hpp>
#include <Creators.hpp>
#include <Print.hpp>

template<class T> inline static int check () {
    Matrix<T> A = randn<T>(3,4), B = A, C = A, D = randn<T>(1);
    T a = D[0];
    std::cout << "A = [" <<std::endl;
    std::cout << A << "];" << std::endl;
    std::cout << a << "*A" <<std::endl;
    std::cout << a*A << std::endl;
    std::cout << "A*" << a  <<std::endl;
    std::cout << A*a << std::endl;
    std::cout << "A*=" << a <<std::endl;
    std::cout << B << std::endl;
    std::cout << "A*A" <<std::endl;
    std::cout << A*A << std::endl;
    std::cout << "A*=A" <<std::endl;
    std::cout << C << std::endl;
    return 0;
}

int main (int args, char** argv) {
    return check<float>() + check<cxfl>() + check<double>() + check<cxdb>();
}
