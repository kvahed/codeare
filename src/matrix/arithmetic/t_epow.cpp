#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

template<class T>
void epow_check () {

	float p = 2.0;
    Matrix<T> A = rand<T>(3,4);
    Matrix<T> C, D;

#ifdef VERBOSE
    std::cout << "A=[\n" << A << "];\n";
    std::cout << "p=" << p << ";\n";
#endif

    C = A^p;
    D = A;
    D ^= p;

#ifdef VERBOSE
    std::cout << "A  = A^p [\n" << C << "];\n";
    std::cout << "A ^= p   [\n" << D << "];\n\n";
#endif

}


int main (int args, char** argv) {
    
    epow_check<float>();
    epow_check<double>();
    epow_check<cxfl>();
    epow_check<cxdb>();
    
    return 0;
    
}
