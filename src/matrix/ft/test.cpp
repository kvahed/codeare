
#include "Algos.hpp"
#include "Creators.hpp"
#include "Params.hpp"

#include <codeare-ft.h>

int main (int args, char** argv) {

    Matrix<cxfl> A = rand<cxfl> (8,8), B;
    DFT<float> ft (8);
    B = ft * A;
    
    return 0;
    
}
