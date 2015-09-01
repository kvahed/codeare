#include <Matrix.hpp>
#include <Creators.hpp>
#include <Print.hpp>

template<class T> int check () {
    Matrix<T> A = randn<T>(4,3), B = A(CR(size(A,0)-1,-1,0),CR()), C = (Matrix<cxdb>) A;
    std::cout << "A = [" << std::endl;
    std::cout << A << std::endl;
    std::cout << "] \n B = A(end:-1:1,:)" << std::endl;
    std::cout << B << std::endl;
    std::cout << "] \n C = (Matrix<T>) A" << std::endl;
    std::cout << C << std::endl;
    return 0;
}

int main (int nargs, const char* argv[]) {
    return check<cxfl>();
}
