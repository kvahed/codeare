#include "Smooth.hpp"

template<class T> int check () {
    Matrix<T> A = randn<T>(10,1);
    A = smooth(A);
    return 0;
}

int main (int, char**) {
    return check<float>();
}
