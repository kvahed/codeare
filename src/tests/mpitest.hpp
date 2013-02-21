#include <math.h>
#include <iostream>
#include "PIO.hpp"

template<> inline
Matrix<float,MPI>::Matrix (const size_t& cols, const size_t& rows) {

    // Get at home
    int info;
    Grid& gd = Workspace::Instance()->GridEnv();
    std::cout << gd;
    
    // Global size
    _bs      = 16;
    _gdim[0] = cols;
    _gdim[1] = rows;
		
    // Local size (only with MPI different from global)
    _dim[0] = numroc_ (&_gdim[0], &_bs, &gd.mc, &izero, &gd.nc);
    _dim[1] = numroc_ (&_gdim[1], &_bs, &gd.mr, &izero, &gd.nr);
	
    // Allocate
    _M.resize(Size());
	
    // Descriptor 
    int dims[2]; dims[0] = _dim[0]; dims[1] = _dim[1];
    descinit_(_desc, &_gdim[0], &_gdim[1], &_bs, &_bs, &izero, &izero, &gd.ct, 
              dims, &info);

//#ifdef BLACS_DEBUG
    printf ("info(%d) desc({%d, %d, %4d, %4d, %d, %d, %d, %d, %4d})\n", 
            info,     _desc[0], _desc[1], _desc[2], _desc[3], 
            _desc[4], _desc[5], _desc[6], _desc[7], _desc[8]);
//#endif
    Cblacs_barrier();
    
}
	
template <class T> bool
mpitest (Connector<T>* rc) {

    Grid& gd = Workspace::Instance()->GridEnv();
    
#ifdef HAVE_MPI

	/* MPI aware matrices */
	Matrix<float,MPI> A (256,256);
	Matrix<float,MPI> B (256,256);

	//A = gd.rk; //operator=(const T& t)
	//B = A;     //operator=(const PMatrix<T>& P)

    //print (B,"B");
    
#else 
	printf ("\nMPI not supported by this build\n");
#endif

	return 0;

}
