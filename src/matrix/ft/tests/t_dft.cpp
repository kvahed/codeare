#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Params.hpp"
#include "DFT.hpp"


int main (int args, char** argv) {

    Matrix<cxfl> A = rand<cxfl> (8,8), B;
    Vector<size_t> sz (2,1);
    sz[0] = 8; sz[1] = 8;
    
    DFT<> ft (sz);
    B = ft * A;
    
    return 0;
    
}
