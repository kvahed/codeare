#include "VXFile.hpp"

using namespace codeare::matrix::io;

template<class T> inline static bool read (const std::string& fname, Matrix<T>& A) {
    VXFile vxf (fname, READ);
    vxf.Read<T>();
    return true;
}

template<class T> inline static bool check () {
  
    Matrix<T> A;
    read("test.dat", A);
    std::cout << wspace << std::endl;
    return true;
  
}

int main (int args, char** argv) {
    if (check<cxfl>())
        return 0;
    else
        return 1;
}
