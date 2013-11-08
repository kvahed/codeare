#ifdef HAVE_SCALAPACK
#include "Matrix.hpp"
#include "Grid.hpp"

template <class T> struct MPITraits {
    Matrix (const size_t cols, const size_t rows) {
        int info;
        Grid& gd = *Grid
        
    }
};

template<> inline
Matrix<float,MPI>::Matrix (const size_t cols, const size_t rows) {

    // Get at home
    int info;
    int izero = 0;
    Grid& gd = *Grid::Instance();
    
    // Global size
    _bs      = 8;
    _gdim[0] = cols;
    _gdim[1] = rows;

    _dim.resize(2);
    _res.resize(2,1.0);
			
    // Local size (only with MPI different from global)
    _dim[0] = numroc_ (&_gdim[0], &_bs, &gd.mc, &izero, &gd.nc);
    _dim[1] = numroc_ (&_gdim[1], &_bs, &gd.mr, &izero, &gd.nr);

    // Allocate
    this->Allocate();
    
    // RAM descriptor 
    int dims[2]; dims[0] = _dim[0]; dims[1] = _dim[1];
    descinit_(_desc, &_gdim[0], &_gdim[1], &_bs, &_bs, &izero,
              &izero, &gd.ct, dims, &info);

#ifdef BLACS_DEBUG
    printf ("info(%d) desc({%d, %d, %4d, %4d, %d, %d, %d, %d, %4d})\n", 
            info,     _desc[0], _desc[1], _desc[2], _desc[3], 
            _desc[4], _desc[5], _desc[6], _desc[7], _desc[8]);
#endif

    char a = 'a';
    Cblacs_barrier (gd.ct, &a);

}

template<> template<> inline Matrix<float,MPI>&
Matrix<float,MPI>::operator= (const Matrix<float,SHM>& M) {
    assert (Size() == M.Size());
    _M = M.Container();
}

template<> template<> inline Matrix<float,MPI>&
Matrix<float,MPI>::operator= (const Matrix<float,MPI>& M) {

    *this = Matrix<float,MPI>(M.Dims());
    _M = M.Container();
    return *this;
    
}

#endif // HAVE_SCALAPACK

