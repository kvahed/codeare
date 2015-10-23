#include <LocalMaxima.hpp>
#include <Trigonometry.hpp>

template<class T> int check () {
    using namespace codeare::matrix::arithmetic;
    Matrix<T> t = linspace<T>(0., 2*PI-1e-5, 512);
    Matrix<T> A = sin(t) + sin(8*t);
    //loc_max_t<T> lm = findLocalMaxima(A);
    Vector<size_t> lm = findLocalMaxima(A);
    return 0;
}

int main (int narg, const char** argv) {
    return check<float>() + check<double>();
}
