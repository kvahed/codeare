#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "IOContext.hpp"
#include "Resample.hpp"
#include "Print.hpp"

template<class T> bool
check_resample () {


	Matrix<T> img = phantom3D<T> (32);
	resample (img,.5,LINEAR);

	return true;

}


int main (int args, char** argv) {

    if (!check_resample<float>())
        return 1;
    if (!check_resample<double>())
    	return 1;
/*    if (!check_interp1<cxfl>())
        return 1;
//    if (!check_ppval<cxdb>())
        return 1;
    if (!check_ppval<short>())
        return 1;
    if (!check_ppval<long>())
    return 1;*/

    return 0;
    
}
