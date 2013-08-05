#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"
#include "Print.hpp"

#define VERBOSE

template<class T, class S> void svd_check () {

    Matrix<T> A = rand<T>(8,3);
    boost::tuple<Matrix<T>,Matrix<S>,Matrix<T> > usv;

#ifdef VERBOSE
    std::cout << "A=\n" << A;
#endif
    
    usv = svd2<T,S> (A);

#ifdef VERBOSE
    std::cout << "U=\n" << boost::get<0>(usv);
    std::cout << "S=\n" << boost::get<1>(usv);
    std::cout << "S=\n" << svd(A);
    std::cout << "V=\n" << boost::get<2>(usv);
    std::cout << std::endl;
#endif

    usv = svd2<T,S> (A, 'S');
#ifdef VERBOSE
    std::cout << "U=\n" << boost::get<0>(usv);
    std::cout << "S=\n" << boost::get<1>(usv);
    std::cout << "V=\n" << boost::get<2>(usv);
    std::cout << std::endl;
#endif

    usv = svd2<T,S> (A, 'A');
#ifdef VERBOSE
    std::cout << "U=\n" << boost::get<0>(usv);
    std::cout << "S=\n" << boost::get<1>(usv);
    std::cout << "V=\n" << boost::get<2>(usv);
    std::cout << std::endl;
#endif

}

int main (int args, char** argv) {

    svd_check<float,float>();
    svd_check<double,double>();
    svd_check<cxfl,float>();
    svd_check<cxdb,double>();
    
    return 0;
    
}
