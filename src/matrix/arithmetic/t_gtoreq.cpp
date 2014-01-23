#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

template<class T>
void emul_check () {

    Matrix<T> A = rand<T>(3,4);
    Matrix<T> B = rand<T>(3,4);
    T a = T(0.0);
    Matrix<cbool> C, D;

#ifdef VERBOSE
    std::cout << "A=[\n" << A << "];\n";
    std::cout << "B=[\n" << B << "];\n";
    std::cout << "a=" << a << ";\n";
#endif

    C = (A >= a);
    D = (a >= A);

#ifdef VERBOSE
    std::cout << "(A >= a) =[\n" << C << "];\n";
    std::cout << "(a >= A) =[\n" << D << "];\n";
    std::cout << std::endl;
#endif

    C = (A >= B);
    D = (B >= A);

#ifdef VERBOSE
    std::cout << "(A >= B) =[\n" << C << "];\n";
    std::cout << "(B >= A) =[\n" << D << "];\n";
    std::cout << std::endl;
#endif
    /*
    D = B;
    B *= a;
    C = B;
    B = D;
    B *= A;

#ifdef VERBOSE
    std::cout << "B*=a[\n" << C << "];\n";
    std::cout << "B*=A[\n" << B << "];\n";
    std::cout << std::endl;
#endif
    */
}


int main (int args, char** argv) {
    
    emul_check<float>();
    emul_check<double>();
    emul_check<cxfl>();
    emul_check<cxdb>();
    
    return 0;
    
}
