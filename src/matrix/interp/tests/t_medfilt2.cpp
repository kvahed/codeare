#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "MedianFilter.hpp"

template<class T> bool
check_medfilt2 () {

    Matrix<T> A = phantom<T>(256) + 2.5e-2 * rand<T>(256);
    A = medfilt2(A,5,5);

    return true;

}


int main (int args, char** argv) {

    if (!check_medfilt2<float>())
        return 1;
    if (!check_medfilt2<double>())
        return 1;
    if (!check_medfilt2<cxfl>())
        return 1;
    if (!check_medfilt2<cxdb>())
        return 1;

    return 0;
    
}
