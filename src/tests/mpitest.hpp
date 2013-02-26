#include "Grid.hpp"
#include "Matrix.hpp"

#include <math.h>
#include <iostream>

template <class T> bool
mpitest (Connector<T>* rc) {

    Grid& gd = *Grid::Instance();
    
#ifdef HAVE_MPI

	/* MPI aware matrices */
	Matrix<float,MPI> A (16,16);
	Matrix<float,MPI> B (16,16);

    //std::cout << A[0] << std::endl;
	//A = 1.0;//gd.rk; //operator=(const T& t)
	//B = A;     //operator=(const PMatrix<T>& P)

#else 
	printf ("\nMPI not supported by this build\n");
#endif

	return 0;

}
