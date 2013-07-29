#include "Grid.hpp"
#include "PMatrix.hpp"
#include "Print.hpp"
#include "PIO.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

#include <math.h>
#include <iostream>



template<class T> bool
check_alloc() {

    Grid& gd = *Grid::Instance();
	Matrix<T,MPI> A (16,16);
    A = rand<T>(size(A));
    return true;

}


int main (int args, char** argv) {

    Grid& gd = *Grid::Instance();
    check_alloc<float>();
    //check_alloc<double>();
    //check_alloc<cxfl>();
    // check_alloc<cxdb>();
    
    return 0;

}
