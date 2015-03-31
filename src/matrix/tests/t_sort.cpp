#include "Matrix.hpp"
#include "Creators.hpp"
#include "Print.hpp"

int main (int narg, char** argv) {
    Matrix<float> A = randn<float>(8,2);
    std::cout << "A" << std::endl;
    std::cout << A << std::endl;
    Vector<size_t> idx = sort (A(Range(),1), DESCENDING);
    std::cout << "idx = sort(A(\":,1\"), DESCENDING)" << std::endl;
    std::cout << idx << std::endl;
    std::cout << "A(\"idx,:\")" << std::endl;
    std::cout << A(idx,Range()) << std::endl;
    idx = sort (A(":,1"));
    std::cout << "idx = sort(A(\":,1\"))" << std::endl;
    std::cout << idx << std::endl;
    std::cout << "A(\"idx,:\")" << std::endl;
    std::cout << A(idx,Range()) << std::endl;
    return 0;
}
