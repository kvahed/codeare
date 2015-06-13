#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"
#include "Print.hpp"

#define VERBOSE

template<class T> void eig_check () {

    Matrix<T> A = rand<T>(4,4);
    typedef typename TypeTraits<T>::CT CT;
    TUPLE <Matrix<T>,Matrix<CT>,Matrix<T> > lev;

#ifdef VERBOSE    
    std::cout << "A =\n" << A << std::endl;
#endif
    
    lev = eig2<T> (A, 'V', 'V');
    
#ifdef VERBOSE    
    std::cout << "ev=\n" << GET<1>(lev) << std::endl;
    std::cout << "ev=\n" << eig(A) << std::endl;
    std::cout << "lv=\n" << GET<0>(lev) << std::endl;
    std::cout << "rv=\n" << GET<2>(lev) << std::endl;
    std::cout << std::endl;
#endif

}

int main (int args, char** argv) {

    eig_check<float>();
    eig_check<double>();
    eig_check<cxfl>();
    eig_check<cxdb>();
    
    return 0;
    
}
