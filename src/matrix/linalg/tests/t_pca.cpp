#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"
#include "Print.hpp"

#define VERBOSE

template<class T> void pca_check () {

    Matrix<T> A = rand<T>(5,3);
    typedef typename TypeTraits<T>::RT RT;
    TUPLE<Matrix<RT>,Matrix<T> > pc;

#ifdef VERBOSE
    std::cout << "A=\n" << A  << std::endl;
#endif
    
    pc = pca2<T> (A);

#ifdef VERBOSE
    std::cout << "C=\n" << GET<0>(pc) << std::endl;
    std::cout << "P=\n" << GET<1>(pc) << std::endl;
    std::cout << std::endl;
#endif

}

int main (int args, char** argv) {

    pca_check<float>();
    pca_check<double>();
    pca_check<cxfl>();
    pca_check<cxdb>();
    
    return 0;
    
}
