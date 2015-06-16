#include "Matrix.hpp"
#include "Creators.hpp"
#include "Print.hpp"

typedef typename Matrix<float>::View::Range R;
typedef typename Matrix<float>::ConstView::Range CR;

int main (int narg, char** argv) {
    Matrix<float> A = randn<float>(8,2);
    std::cout << "A" << std::endl;
    std::cout << A << std::endl;
    Vector<size_t> idx = sort (Matrix<float>(A(CR(),CR(1,1))), DESCENDING);
    std::cout << "idx = sort(A(\":,1\"), DESCENDING)" << std::endl;
    std::cout << idx << std::endl;
    std::cout << "A(\"idx,:\")" << std::endl;
    std::cout << Matrix<float>(A(CR(idx),CR())) << std::endl;
    return 0;
}
