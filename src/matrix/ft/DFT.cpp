#include "DFT.hpp"

//template Matrix<cxfl> fftshift (const Matrix<cxfl>&, const bool);
//template Matrix<cxdb> ifftshift (const Matrix<cxdb>&);

template class DFT<float>;
template class DFT<double>;
