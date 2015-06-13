#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"
#include "Print.hpp"

#define VERBOSE

template<class T> void svd_check () {

    Matrix<T> A = rand<T>(5,3);
    typedef typename TypeTraits<T>::RT RT;
    TUPLE<Matrix<T>,Matrix<RT>,Matrix<T> > usv;

#ifdef VERBOSE
    std::cout << "A=\n" << A  << std::endl;
#endif
    
    usv = svd2<T> (A);

#ifdef VERBOSE
    std::cout << "U=\n" << GET<0>(usv) << std::endl;
    std::cout << "S=\n" << GET<1>(usv) << std::endl;
    std::cout << "S=\n" << svd(A) << std::endl;
    std::cout << "V=\n" << GET<2>(usv) << std::endl;
    std::cout << std::endl;
#endif

    usv = svd2<T> (A, 'S');
#ifdef VERBOSE
    std::cout << "U=\n" << GET<0>(usv) << std::endl;
    std::cout << "S=\n" << GET<1>(usv) << std::endl;
    std::cout << "V=\n" << GET<2>(usv) << std::endl;
    std::cout << std::endl;
#endif

    usv = svd2<T> (A, 'A');
#ifdef VERBOSE
    std::cout << "U=\n" << GET<0>(usv) << std::endl;
    std::cout << "S=\n" << GET<1>(usv) << std::endl;
    std::cout << "V=\n" << GET<2>(usv) << std::endl;
    std::cout << std::endl;
#endif

}

int main (int args, char** argv) {

    svd_check<float>();
    svd_check<double>();
    svd_check<cxfl>();
    svd_check<cxdb>();
    
    return 0;
    
}
