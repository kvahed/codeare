#include <Matrix.hpp>
#include <Creators.hpp>
#include <Print.hpp>

template<class T> inline static int check () {
    Matrix<T> A = randn<T>(4,3);
    std::cout << "A = [" << std::endl;
    std::cout << A << std::endl;
    std::cout << "];"<<std::endl;
    std::cout << "sum(A) " << std::endl;
    std::cout << sum(A) << std::endl;
    std::cout << "sum(A,1) " << std::endl;
    std::cout << sum(A,1) << std::endl << std::endl;
    return 0;
}

int main (int args, char** argv) {
    return check<float>() + check<double>() + check<cxfl>() + check<cxdb>();
}
