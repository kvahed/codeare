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
	dot(x,y), b = dotc(x,y), c = 0., d = 0.; 

	for (size_t i = 0; i < n; ++i) {
		c += x[i]*y[i];
		d += TypeTraits<T>::Conj(x[i])*y[i];
	}

#ifdef VERBOSE
    std::cout << "x=\n" << x << "y=\n" << y;
    std::cout << "y * y=\n" << a << " " << c << std::endl;
    std::cout << "x**H * y=\n" << b << " " << d << "\n" << std::endl;
#endif
	return (TypeTraits<T>::Abs((a-c)/a)<1.e-6 && 
			TypeTraits<T>::Abs((b-d)/b)<1.e-6);
}

int main (int args, char** argv) {
	return (dot_check<cxfl>() && dot_check<cxdb>() && 
			dot_check<float>() && dot_check<double>()) ? 0 : 1;
}
