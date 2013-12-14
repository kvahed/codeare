#include "MRI.hpp"

template Matrix<float> phase_combine (const Matrix<float>& A, const size_t dim);
template Matrix<double> phase_combine (const Matrix<double>& A, const size_t dim);
template Matrix<cxfl> phase_combine (const Matrix<cxfl>& A, const size_t dim);
template Matrix<cxdb> phase_combine (const Matrix<cxdb>& A, const size_t dim);
