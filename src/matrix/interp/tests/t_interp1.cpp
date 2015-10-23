#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Interpolate.hpp"
#include "Print.hpp"

//#define VERBOSE

template<class T> bool
check_interp1 () {

	size_t n = 8;
	size_t m = 5;

	Matrix<double> x(n+1,1),xi (n*m+1,1);
	Matrix<T> y, yi;
	y = randn<T>(n+1,1);

	for (size_t i = 0; i <= n; ++i)
		x[i] = (double)i;
	for (size_t i = 0; i <= n*m; ++i) 
		xi[i] = (double)i*(1./m);

	yi = interp1(x, y, xi, INTERP::AKIMA);

	std::cout << yi << "\n";

	return true;

}


int main (int args, char** argv) {

    if (!check_interp1<float>())
        return 1;
    if (!check_interp1<double>())
    	return 1;
    if (!check_interp1<cxfl>())
        return 1;
    if (!check_interp1<cxdb>())
        return 1;

    return 0;
    
}
