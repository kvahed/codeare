#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "IOContext.hpp"
#include "MedianFilter.hpp"

template<class T> bool
check_medfilt2 () {

    IOContext f = fopen ("test.mat",WRITE);

    Matrix<T> A = phantom<T>(256) + 2.5e-2 * rand<T>(256);
    fwrite (f, A, "A");
    A = medfilt2(A);
    fwrite (f, A, "A_filt");

    fclose (f);

    return true;

}


int main (int args, char** argv) {

    if (!check_medfilt2<cxfl>())
        return 1;

    return 0;
    
}
