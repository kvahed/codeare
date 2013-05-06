#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Params.hpp"
#include "DWT.hpp"


int main (int args, char** argv) {

    size_t sl = 32;
    
    Matrix<cxfl> A = rand<cxfl> (32,32), B;
    
    DWT<cxfl> wt (sl);
    B = wt * A;
    
    return 0;
    
}
