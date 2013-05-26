#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"
#include "Print.hpp"

template<class T, class S> void eig_check () {

    Matrix<T> A = rand<T>(6,6);
    boost::tuple <Matrix<T>,Matrix<S>,Matrix<T> > lev;

#ifdef VERBOSE    
    std::cout << "A =\n" << A;
#endif
    
    lev = eig2<T,S> (A, 'V', 'V');
    
#ifdef VERBOSE    
    std::cout << "ev=\n" << boost::get<1>(lev);
    std::cout << "ev=\n" << eig(A);
    std::cout << "lv=\n" << boost::get<0>(lev);
    std::cout << "rv=\n" << boost::get<2>(lev);
    std::cout << std::endl;
#endif

}

int main (int args, char** argv) {

    eig_check<float,cxfl>();
    eig_check<double,cxdb>();
    eig_check<cxfl,cxfl>();
    eig_check<cxdb,cxdb>();
    
    return 0;
    
}
