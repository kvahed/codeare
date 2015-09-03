#include <Matrix.hpp>
#include <Creators.hpp>
#include <Algos.hpp>
#include <Print.hpp>

template<class T> int check () {
    Matrix<T> A = randn<T>(2,3), B = flipud(A), C = fliplr(A);
    std::cout << "A = [" << std::endl;
    std::cout << A << std::endl;
    std::cout << "]\n flipud(A)" << std::endl;
    std::cout << flipud(A) << std::endl;
    std::cout << "]\n fliplr(A)" << std::endl;
    std::cout << fliplr(A) << std::endl;
    return 0;
}

int main (int args, char** argv) {

    return check<float>() + check<double>() + check<cxfl>() + check<cxdb>();
    
}

