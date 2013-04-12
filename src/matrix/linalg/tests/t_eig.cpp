#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"

template<class T, class S> void eig_check () {

    Matrix<T> A = rand<T>(6,6), lv(6,6), rv(6,6);
    Matrix<S> ev;

#ifdef VERBOSE    
    std::cout << "A =\n" << A;
#endif
    
    eig (A, ev, lv, rv, 'N', 'N');
    
#ifdef VERBOSE    
    std::cout << "ev=\n" << ev;
    std::cout << "lv=\n" << lv;
    std::cout << "rv=\n" << rv << std::endl;
#endif

}

int main (int args, char** argv) {

    //eig_check<float,cxfl>();
    //eig_check<double,cxdb>();
    //eig_check<cxfl,cxfl>();
    //eig_check<cxdb,cxdb>();
    
    return 0;
    
}
