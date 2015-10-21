#include <Smooth.hpp>
#include <Print.hpp>
#include <Trigonometry.hpp>

using namespace codeare::matrix::arithmetic;

template<class T> int check () {
	typedef typename TypeTraits<T>::RT RT;
	Matrix<RT> t = linspace<T>(0,2*PI,100);
	t = squeeze(5*sin(t)+1*sin(8*t));
    Matrix<T> A = randn<T>(100,1)+t;
    std::cout << "A = [" <<std::endl;
    std::cout << A << "];" << std::endl;
    A = smooth(A, 9);
    std::cout << "smooth(A,9)" <<std::endl;
    std::cout << A << std::endl;
    return 0;
}

int main (int, char**) {
    return check<double>() + check<float>();
}
