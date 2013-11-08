#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Params.hpp"
#include "DFT.hpp"


int main (int args, char** argv) {

    Matrix<cxfl> A = phantom<cxfl> (32), B;
    Matrix<size_t> sz (2,1);
    sz[0] = 32; sz[1] = 32;
    
    DFT<float> ft (sz);
    B = ft * A;

    std::cout << A << std::endl;
    std::cout << B << std::endl;
    
    return 0;
    
}
