#include <Matrix.hpp>
#include <Creators.hpp>
#include <Print.hpp>

template<class S, class T> int check () {
    Matrix<T> A = randn<T>(3,4);
    Matrix<S> B = randn<S>(1,1);
    S s = B[0];
    std::cout << "A = [" << std::endl;
    std::cout << A << std::endl;
    std::cout << "]\n"<< s<< "+A" << std::endl;
    std::cout << s + A << std::endl;
    std::cout << s << "*A" << std::endl;
    std::cout << s * A << std::endl;
    std::cout << s << "./A" << std::endl;
    std::cout << s / A << std::endl;
    std::cout << s << "-A" << std::endl;
    std::cout << s - A << std::endl;
    return 0;
}

int main (int args, const char* argv[]) {
    return check<float,cxfl>() + check<double,float>();
}
