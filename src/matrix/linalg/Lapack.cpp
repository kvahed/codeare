#include "Lapack.hpp"

template Matrix<float> chol (const Matrix<float>& A, const char& = 'U');
template Matrix<double> chol (const Matrix<double>& A, const char& = 'U');
template Matrix<cxfl> chol (const Matrix<cxfl>& A, const char& = 'U');
template Matrix<cxdb> chol (const Matrix<cxdb>& A, const char& = 'U');

template int eig (const Matrix<cxfl>&, Matrix<cxfl>&, Matrix<cxfl>&, Matrix<cxfl>&, const char& = 'N', const char& = 'N');
template int eig (const Matrix<cxdb>&, Matrix<cxdb>&, Matrix<cxdb>&, Matrix<cxdb>&, const char& = 'N', const char& = 'N');
template int eig (const Matrix<float>&, Matrix<float>&, Matrix<float>&, Matrix<float>&, const char& = 'N', const char& = 'N');
template int eig (const Matrix<double>&, Matrix<double>&, Matrix<double>&, Matrix<double>&, const char& = 'N', const char& = 'N');

template float dot (const Matrix<float>&, const Matrix<float>&);
template double dot (const Matrix<double>&, const Matrix<double>&);
template cxfl dot (const Matrix<cxfl>&, const Matrix<cxfl>&);
template cxdb dot (const Matrix<cxdb>&, const Matrix<cxdb>&);

template float dotc (const Matrix<float>&, const Matrix<float>&);
template double dotc (const Matrix<double>&, const Matrix<double>&);
template cxfl dotc (const Matrix<cxfl>&, const Matrix<cxfl>&);
template cxdb dotc (const Matrix<cxdb>&, const Matrix<cxdb>&);

template Matrix<float> gemm (const Matrix<float>&, const Matrix<float>&, const char& = 'N', const char& = 'N');
template Matrix<double> gemm (const Matrix<double>&, const Matrix<double>&, const char& = 'N', const char& = 'N');
template Matrix<cxfl> gemm (const Matrix<cxfl>&, const Matrix<cxfl>&, const char& = 'N', const char& = 'N');
template Matrix<cxdb> gemm (const Matrix<cxdb>&, const Matrix<cxdb>&, const char& = 'N', const char& = 'N');

template Matrix<float> gemv (const Matrix<float>&, const Matrix<float>&, const char& = 'N');
template Matrix<double> gemv (const Matrix<double>&, const Matrix<double>&, const char& = 'N');
template Matrix<cxfl> gemv (const Matrix<cxfl>&, const Matrix<cxfl>&, const char& = 'N');
template Matrix<cxdb> gemv (const Matrix<cxdb>&, const Matrix<cxdb>&, const char& = 'N');

template Matrix<float> inv (const Matrix<float>&);
template Matrix<double> inv (const Matrix<double>&);
template Matrix<cxfl> inv (const Matrix<cxfl>&);
template Matrix<cxdb> inv (const Matrix<cxdb>&);

template double norm (const Matrix<float>&);
template double norm (const Matrix<double>&);
template double norm (const Matrix<cxfl>&);
template double norm (const Matrix<cxdb>&);

template Matrix<float> pinv (const Matrix<float>&, const char& = 'N');
template Matrix<double> pinv (const Matrix<double>&, const char& = 'N');
template Matrix<cxfl> pinv (const Matrix<cxfl>&, const char& = 'N');
template Matrix<cxdb> pinv (const Matrix<cxdb>&, const char& = 'N');

template int svd (const Matrix<float>&, Matrix<float>&, Matrix<float>&, Matrix<float>&, const char& = 'N');
template int svd (const Matrix<double>&, Matrix<double>&, Matrix<double>&, Matrix<double>&, const char& = 'N');
template int svd (const Matrix<cxfl>&, Matrix<float>&, Matrix<cxfl>&, Matrix<cxfl>&, const char& = 'N');
template int svd (const Matrix<cxdb>&, Matrix<double>&, Matrix<cxdb>&, Matrix<cxdb>&, const char& = 'N');

