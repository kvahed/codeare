#include <Matrix.hpp>
#include <Creators.hpp>
#include <Print.hpp>

template<class T> inline static int check () {
    Matrix<T> A (1,2,1,3,1,4,1,5,1,6,1);
    std::cout << "size(A) = " << std::endl;
    std::cout << size(A) << std::endl;
    std::cout << "size(squeeze(A)) = " << std::endl;
    std::cout << size(squeeze(A)) << std::endl;
    return 0;
}

int main (int args, char** argv) {
    return check<float>();
}
