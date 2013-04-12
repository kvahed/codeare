#include "DFT.hpp"

template Matrix<cxfl> fftshift (const Matrix<cxfl>&);
template Matrix<cxdb> fftshift (const Matrix<cxdb>&);

template class DFT<float>;
template class DFT<double>;
