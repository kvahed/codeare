#include "Lapack.hpp"
#include "Matrix.hpp"

template class Matrix<float,SHM>;
template class Matrix<double,SHM>;
template class Matrix<cxfl,SHM>;
template class Matrix<cxdb,SHM>;
template class Matrix<bool,SHM>;

