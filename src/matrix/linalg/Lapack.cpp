#include "Lapack.hpp"

template Matrix<float> chol (const Matrix<float>& A, const char& = 'U');
template Matrix<double> chol (const Matrix<double>& A, const char& = 'U');
template Matrix<cxfl> chol (const Matrix<cxfl>& A, const char& = 'U');
template Matrix<cxdb> chol (const Matrix<cxdb>& A, const char& = 'U');

template boost::tuple<Matrix<float>,Matrix<cxfl>,Matrix<float> > eig2<float,cxfl> (const Matrix<float>&, const char& = 'N', const char& = 'N');
template boost::tuple<Matrix<double>,Matrix<cxdb>,Matrix<double> > eig2<double,cxdb> (const Matrix<double>&, const char& = 'N', const char& = 'N');
template boost::tuple<Matrix<cxfl>,Matrix<cxfl>,Matrix<cxfl> > eig2<cxfl,cxfl> (const Matrix<cxfl>&, const char& = 'N', const char& = 'N');
template boost::tuple<Matrix<cxdb>,Matrix<cxdb>,Matrix<cxdb> > eig2<cxdb,cxdb> (const Matrix<cxdb>&, const char& = 'N', const char& = 'N');

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

template boost::tuple<Matrix<float>,Matrix<float>,Matrix<float> > svd2<float,float> (const Matrix<float>&, const char& = 'N');
template boost::tuple<Matrix<double>,Matrix<double>,Matrix<double> > svd2<double,double> (const Matrix<double>&, const char& = 'N');
template boost::tuple<Matrix<cxfl>,Matrix<float>,Matrix<cxfl> > svd2<cxfl,float> (const Matrix<cxfl>&, const char& = 'N');
template boost::tuple<Matrix<cxdb>,Matrix<double>,Matrix<cxdb> > svd2<cxdb,double> (const Matrix<cxdb>&, const char& = 'N');

