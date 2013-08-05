#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "IOContext.hpp"
#include "PolyVal.hpp"

#define VERBOSE

template<class T> bool
check_ppval () {

	size_t n = 3;
	size_t m = 20;

	Matrix<double> x (n+1,1);
	Matrix<T> y = randn<T>(n+1,1), yi (n*m+1,1);

	for (size_t i = 0; i <= n; ++i)
		x[i] = (double)i;

	PolyVal<T> pv = PolyVal<T>(x, y, INTERP::CSPLINE);
	
	for (size_t i = 0; i <= n*m; ++i)
		yi[i] = pv.Lookup((double)i*(1./m));
#ifdef VERBOSE
		std::cout << yi << "\n";
#endif
	return true;

}


int main (int args, char** argv) {

    if (!check_ppval<float>())
        return 1;
    if (!check_ppval<double>())
    	return 1;
    if (!check_ppval<cxfl>())
        return 1;
    if (!check_ppval<cxdb>())
        return 1;

    return 0;
    
}
