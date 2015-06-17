#include "Matrix.hpp"
#include "Creators.hpp"
#include "Print.hpp"

int main (int narg, char** argv) {
    Matrix<float> A = randn<float>(8,1);
    std::cout << "A" << std::endl;
    std::cout << A << std::endl;
    Vector<size_t> idx = sort (A, DESCENDING);
    std::cout << "idx = sort(A(\":,1\"), DESCENDING)" << std::endl;
    std::cout << idx << std::endl;
    return 0;
}
