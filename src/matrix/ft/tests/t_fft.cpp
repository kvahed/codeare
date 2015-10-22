#include "DFT.hpp"
#include "IOContext.hpp"

template<class T> inline int check () {
    Matrix<cxfl> A = phantom<cxfl>(9);
    std::cout << "A = [" <<std::endl;
    std::cout << A << "];" << std::endl;
    std::cout << "fft(A)"  <<std::endl;
    std::cout << fft(A) << std::endl;
    std::cout << "fft(A,0,false)"  <<std::endl;
    std::cout << fft(A,0,false) << std::endl;
    std::cout << "fft(A,1,false)"  <<std::endl;
    std::cout << fft(A,1,false) << std::endl;
    return 0;
}

int main (int narg, char** argv) {
    return check<cxfl>() + check<cxdb>();
    return 0;
}
