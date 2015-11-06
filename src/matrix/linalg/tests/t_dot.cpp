#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"
#include "Print.hpp"

#define VERBOSE

template<class T> bool dot_check () {

	size_t n = 4;
    Matrix<T> x = rand<T> (n,1);
    Matrix<T> y = rand<T> (n,1);
#ifdef VERBOSE
    std::cout << "x=\n" << x  << std::endl;
    std::cout << "y=\n" << y  << std::endl;
#endif
	T a = dot(x,y), c = 0.;
	for (size_t i = 0; i < n; ++i)
		c += x[i]*y[i];
#ifdef VERBOSE
    std::cout << "y * y=\n" << a << " " << c << std::endl;
#endif
    T b = dotc(x,y), d = 0.; 

	for (size_t i = 0; i < n; ++i) 
		d += TypeTraits<T>::Conj(x[i])*y[i];
#ifdef VERBOSE
    std::cout << "x**H * y=\n" << b << " " << d << "\n" << std::endl;
#endif
    std::cout << TypeTraits<T>::Abs((a-c)/a) << " " << TypeTraits<T>::Abs((b-d)/b) << std::endl;
	return (TypeTraits<T>::Abs((a-c)/a)<1.e-4 && 
			TypeTraits<T>::Abs((b-d)/b)<1.e-4);
}

int main (int args, char** argv) {
	return (dot_check<float>() && dot_check<double>() && 
            dot_check<cxfl>() && dot_check<cxdb>()) ? 0 : 1;
}
