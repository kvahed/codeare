#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Params.hpp"
#include "DFT.hpp"


int main (int args, char** argv) {

    Matrix<cxfl> A = rand<cxfl> (8,8), B;

    DFT<cxfl> ft (size(A));
    B = ft * A;
    
    return 0;
    
}
