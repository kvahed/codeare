/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but 
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301  USA
 */

#ifndef __LAPACK_HPP__
#define __LAPACK_HPP__

#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "LapackTraits.hpp"
#include "CX.hpp"

#include <boost/tuple/tuple.hpp>

/**
 * @brief         Eigenvalue decomposition
 * 
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl> m = rand<cxfl> (10,10), ev;
 *   Matrix<float> lv, rv;
 *   boost::tuple<Matrix<float>,Matrix<cxfl>,Matrix<float>> ler
 *   ler = eig2 (m);
 *   lv = boost::get<0>(ler)
 *   ev = boost::get<1>(ler);
 *   rv = boost::get<2>(ler)
 * @endcode
 * where m is the complex decomposed matrix, lv and rv are the left hand and right 
 * hand side eigen vectors and ev is the eigenvalue vector.
 *
 * @see           LAPACK driver xGEEV
 *
 * @param  m      Matrix for decomposition
 * @param  jobvl  Compute left vectors ('N'/'V')
 * @param  jobvr  Compute right vectors ('N'/'V')
 * @return        Eigenvectors and values
 */
template <class T, class S> inline boost::tuple< Matrix<T>, Matrix<S>, Matrix<T> >
eig2 (const Matrix<T>& m, char jobvl = 'N', char jobvr = 'N') {
    
    typedef typename LapackTraits<T>::RType T2;
    T t = (T)0;
    
    boost::tuple< Matrix<T>, Matrix<S>, Matrix<T> > ret;
    
    
    assert (jobvl == 'N' || jobvl =='V');
    assert (jobvr == 'N' || jobvr =='V');
    assert (issquare(m));
    
    int    N     =  size(m, 0);
    int    lda   =  N;
    int    ldvl  = (jobvl == 'V') ? N : 1;
    int    ldvr  = (jobvr == 'V') ? N : 1;
    int    info  =  0;
    int    lwork = -1;
    
    // Appropriately resize the output
    boost::get<1>(ret) = Matrix<S>(N,1);
    
    if (jobvl == 'V') boost::get<0>(ret) = Matrix<T>(N,N);
    if (jobvr == 'V') boost::get<2>(ret) = Matrix<T>(N,N);
    
    Matrix<S>& ev = boost::get<1>(ret);
    Matrix<T> &lv = boost::get<0>(ret), &rv = boost::get<2>(ret);


    // Workspace for real numbers
    container<T2> rwork = (is_complex(t)) ?
        container<T2>(2*N) : container<T2>(1);

    // Workspace for complex numbers
    container<T>  work = container<T>(1);
    
    // Need copy. Lapack destroys A on output.
    Matrix<T> a = m;
    
    // Work space query
    LapackTraits<T>::geev (jobvl, jobvr, N, &a[0], lda, &ev[0], &lv[0], ldvl, &rv[0], ldvr, &work[0], lwork, &rwork[0], info);
    
    // Initialize work space
    lwork = (int) TypeTraits<T>::Real (work[0]);
    work.resize(lwork);
    
    // Actual Eigenvalue computation
    LapackTraits<T>::geev (jobvl, jobvr, N, &a[0], lda, &ev[0], &lv[0], ldvl, &rv[0], ldvr, &work[0], lwork, &rwork[0], info);

    if (info > 0) {
        printf ("\nERROR - XGEEV: the QR algorithm failed to compute all the\n eigenvalues, and no eigenvectors have " \
                "been computed;\n elements %d+1:N of ev contain eigenvalues which\n have converged.\n\n", info) ;
    } else if (info < 0)
        printf ("\nERROR - XGEEV: the %d-th argument had an illegal value.\n\n", -info);
    
    return ret;
    
}

inline Matrix<cxfl> eig (const Matrix<float>& m) {
    return boost::get<1>(eig2<float,cxfl>(m,'N','N'));
}
inline Matrix<cxdb> eig (const Matrix<double>& m) {
    return boost::get<1>(eig2<double,cxdb>(m,'N','N'));
}
inline Matrix<cxfl> eig (const Matrix<cxfl>& m) {
    return boost::get<1>(eig2<cxfl,cxfl>(m,'N','N'));
}
inline Matrix<cxdb> eig (const Matrix<cxdb>& m) {
    return boost::get<1>(eig2<cxdb,cxdb>(m,'N','N'));
}


/**
 * @brief           Singular value decomposition
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl>  m = rand<cxfl> (20,10), u, v;
 *   Matrix<float> s;
 *   boost::tuple<Matrix<cxfl>, Matrix<float>, Matrix<cxfl>> usv;
 *   usv = svd2 (m);
 *   u = boost::get<0>(usv);
 *   s = boost::get<1>(usv);
 *   v = boost::get<2>(usv);
 * @endcode
 * where m is the complex decomposed matrix, u and v are the left hand and right 
 * hand side singular vectors and s is the vector of real singular values.
 *
 * @see             LAPACK driver xGESDD
 * 
 * @param  m       Incoming matrix
 * @param  jobz     Computation mode<br/>
 *                  'A': all m columns of U and all n rows of VT are returned in the arrays u and vt<br/>
 *                  'S', the first min(m, n) columns of U and the first min(m, n) rows of VT are returned in the arrays
 *                       u and vt;<br/>
 *                  'O', then<br/>
 *                  &nbsp;&nbsp;&nbsp;&nbsp;if m >= n, the first n columns of U are overwritten in the array a and all
 *                       rows of VT are returned in the array vt;<br/>
 *                  &nbsp;&nbsp;&nbsp;&nbsp;if m < n, all columns of U are returned in the array u and the first m rows
 *                       of VT are overwritten in the array a;<br/>
 *                  'N', no columns of U or rows of VT are computed (default).
 * @return          Signular vectors and values
 */

template<class T, class S> inline boost::tuple<Matrix<T>,Matrix<S>,Matrix<T> >
svd2 (const Matrix<T>& M, char jobz = 'N') {
    
    typedef typename LapackTraits<T>::RType T2;
    boost::tuple<Matrix<T>,Matrix<S>,Matrix<T> > ret;
    
    assert (is2d(M));
    assert (jobz == 'N' || jobz == 'S' || jobz == 'A');
    
    Matrix<T> A (M);
    
    int   m, n, lwork, info = 0, lda, mn, ldu = 1, ucol = 1, ldvt = 1, vtcol = 1;

    container<T2> rwork;
    container<T>  work = container<T>(1);
    
    m     =  A.Height(); n = A.Width();
    lwork = -1;
    lda   =  m;
    mn    =  MIN(m,n);
    
    if (jobz != 'N')
        ldu = m;
    
    if      (jobz == 'A') {
        ldvt =  n;
        vtcol = (m>=n) ? mn : n;
        ucol =  m;
    } else if (jobz == 'S') {
        ucol = mn;
        ldvt = mn;
        vtcol = (m>=n) ? mn : n;
    }

    Matrix<S>& s = boost::get<1>(ret) = Matrix<S>(  mn,  IONE);
    Matrix<T>& U = boost::get<0>(ret) = Matrix<T>( ldu,  ucol);
    Matrix<T>& V = boost::get<2>(ret) = Matrix<T>(ldvt, vtcol);
    
    container<int> iwork = container<int>(8 * mn);
    
    size_t nr = (typeid(T) == typeid(cxfl) || typeid(T) == typeid(cxdb)) ?
            ((jobz == 'N') ? mn * 7 : mn * (5 * mn + 7)) : 1;
    rwork.resize(nr);

    // Workspace query
    LapackTraits<T>::gesdd (jobz, m, n, &A[0], lda, &s[0], &U[0], ldu, &V[0],
            ldvt, &work[0], lwork, &rwork[0], &iwork[0], info);
    
    lwork = (int) TypeTraits<T>::Real(work[0]);
    work.resize(lwork);
    
    // SVD
    LapackTraits<T>::gesdd (jobz, m, n, &A[0], lda, &s[0], &U[0], ldu, &V[0],
            ldvt, &work[0], lwork, &rwork[0], &iwork[0], info);
    
    
    // Traspose the baby
    V = !V;
    V = conj(V);
    
    if (info > 0)
        printf ("\nERROR - XGESDD: The updating process of SBDSDC did not converge.\n\n");
    else if (info < 0)
        printf ("\nERROR - XGESDD: The %i-th argument had an illegal value.\n\n", -info); 
    
    return ret;
    
} 
// Convenience calls (for s = svd (A))    
inline Matrix<float> svd (const Matrix<float>& A) {
    return boost::get<1>(svd2<float,float>(A));
}
inline Matrix<double> svd (const Matrix<double>& A) {
    return boost::get<1>(svd2<double,double>(A));
}
inline Matrix<float> svd (const Matrix<cxfl>& A) {
    return boost::get<1>(svd2<cxfl,float>(A));
} 
inline Matrix<double> svd (const Matrix<cxdb>& A) {
    return boost::get<1>(svd2<cxdb,double>(A));
}

/**
 * @brief                Invert quadratic well conditioned matrix
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl> m = rand<cxfl> (10,10);
 *  
 *   m = inv (m);
 * @endcode
 *
 * @see                  Lapack xGETRF/xGETRI
 * 
 * @param  m             Matrix
 * @return               Inverse
 */
template <class T> inline Matrix<T> 
inv (const Matrix<T>& m) {
    
    // 2D 
    assert(issquare(m));
    
    int N = (int) size (m,0);    
    Matrix<T> res = m;

    int  info = 0;
    container<int> ipiv = container<int>(N);
    
    // LU Factorisation -------------------
    LapackTraits<T>::getrf (N, N, &res[0], N, &ipiv[0], info);
    
    if (info < 0)
        printf ("\nERROR - DPOTRI: the %i-th argument had an illegal value.\n\n", -info);
    else if (info > 1)
        printf ("\nERROR - DPOTRI: the (%i,%i) element of the factor U or L is\n zero, and the inverse could not be " \
                "computed.\n\n", info, info);
    
    int lwork = -1; 
    container<T> work = container<T>(1);
    
    // Workspace determination ------------
    LapackTraits<T>::getri (N, &res[0], N, &ipiv[0], &work[0], lwork, info);
    
    // Work memory allocation -------------
    lwork = (int) TypeTraits<T>::Real (work[0]);
    work.resize(lwork);
    
    // Inversion --------------------------
    LapackTraits<T>::getri (N, &res[0], N, &ipiv[0], &work[0], lwork, info);

    if (info < 0)
        printf ("\nERROR - XGETRI: The %i-th argument had an illegal value.\n\n", -info);
    else if (info > 0)
        printf ("\nERROR - XGETRI: The leading minor of order %i is not\n positive definite, and the factorization could " \
                "not be\n completed.", info);
    
    return res;
    
} 
    

/**
 * @brief                Moore penrose pseudo-invert
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl> m = rand<cxfl> (10,6);
 *
 *   m = pinv (m);
 * @endcode
 *
 * @see                  Lapack xGELS
 *
 * @param  m             Matrix
 * @param  trans         Transpose m before pinv?
 * @return               Pseudo-inverse
 */
template<class T> inline Matrix<T>
pinv (const Matrix<T>& m, char trans = 'N') {
    
    Matrix<T> mm (m);
    
    assert (is2d(m));

    container<T> work = container<T>(1);
    
    int  M      =  size(m, 0);
    int  N      =  size(m, 1);
    
    int  nrhs   =  M;
    int  lda    =  M;
    int  ldb    =  MAX(M,N);
    int  lwork  = -1;
    int  info   =  0;

    Matrix<T> b = eye<T>(ldb);
    
    LapackTraits<T>::gels (trans, M, N, nrhs, mm, lda, b, ldb, work, lwork, info);
    
    lwork = (int) TypeTraits<T>::Real(work[0]);
    work.resize(lwork);

    LapackTraits<T>::gels (trans, M, N, nrhs, mm, lda, b, ldb, work, lwork, info);
    
    if (M > N)
        for (int i = 0; i < M; i++)
            memcpy (&b[i*N], &b[i*M], N * sizeof(T));
    
    b = resize (b, N, M);

    if (info > 0)
        printf ("ERROR XGELS: the algorithm for computing the SVD failed to converge;\n %i off-diagonal elements " \
                "of an intermediate bidiagonal form\n did not converge to zero.", info);
    else if (info < 0)
        printf ("ERROR XGELS: the %i-th argument had an illegal value.", -info);
    
    return b;
    
}


/**
 * @brief        Cholesky decomposition of positive definite quadratic matrix
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl> m = rand<cxfl> (20,10);
 *   
 *   m = m.prodt(m); // m*m' Must be positive definite 
 *   m = chol (m);
 * @endcode
 *
 * @see          LAPACK driver xPOTRF
 * 
 * @param  A     Incoming matrix
 * @param  uplo  Use upper/lower triangle for decomposition ('U': default/'L')
 * @return       Cholesky decomposition
 */
template<class T> inline Matrix<T> 
chol (const Matrix<T>& A, char uplo = 'U') {
    
    assert(is2d(A));
    
    Matrix<T> res  = A;
    int       info = 0, n = A.Height();
    
    LapackTraits<T>::potrf (uplo, n, &res[0], n, info);
    
    if (info > 0)
        printf ("\nERROR - XPOTRF: the leading minor of order %i is not\n positive definite, and the factorization " \
                "could not be\n completed!\n\n", info);
    else if (info < 0)
        printf ("\nERROR - XPOTRF: the %i-th argument had an illegal value.\n\n!", -info);
    
    for (size_t i = 0; i < n-1; i++)
        for (size_t j = i+1; j < n; j++)
            res(j,i) = T(0);
    
    return res;
    
}


/**
 * @brief          Matrix matrix multiplication
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl> m = rand<cxfl> (20,10);
 *   Matrix<cxfl> x = rand<cxfl> (10, 6);
 *  
 *   m   = gemm (m, b, 'N', 'C');
 * @endcode
 *
 * @see            BLAS routine xGEMM
 *
 * @param  A       Left factor
 * @param  B       Right factor
 * @param  transa  (N: A*... | T: A.'*... | C: A'*...) transpose left factor
 * @param  transb  (N: ...*B | T: ...*B.' | C: ...*B') transpose right factor
 * @return         Product
 */
template<class T> inline Matrix<T> 
gemm (const Matrix<T>& A, const Matrix<T>& B, char transa = 'N', char transb = 'N') {
    
    assert (isvec(A)||is2d(A));
    assert (isvec(B)||is2d(B));
    
	int aw, ah, bw, bh, m, n, k;
	T   alpha = (T)1., beta = (T)0.;
    
    aw = (int)size(A,1); ah = (int)size(A,0), bw = (int)size(B,1), bh = (int)size(B,0);
    
    // Check inner dimensions
	if      ( transa == 'N'                   &&  transb == 'N'                  ) assert (aw == bh);
	else if ( transa == 'N'                   && (transb == 'T' || transb == 'C')) assert (aw == bw);
	else if ((transa == 'T' || transa == 'C') &&  transb == 'N'                  ) assert (ah == bh);
	else if ((transa == 'T' || transa == 'C') && (transb == 'T' || transb == 'C')) assert (ah == bw);
	
	if (transa == 'N') {
		m = ah;
		k = aw;
	} else if (transa == 'T' || transa == 'C') {
		m = aw;
		k = ah;
	}
	
	if (transb == 'N')
		n = bw;
	else if (transb == 'T' || transb == 'C')
		n = bh;

	Matrix<T> C(m,n);
	
	LapackTraits<T>::gemm (transa, transb, m, n, k, alpha, A.Ptr(), ah, B.Ptr(), bh, beta, &C[0], m);
    
	return C;
	
}

template<class T, paradigm P, const bool& b> Matrix<T,P>
Matrix<T,P,b>::prod  (const Matrix<T,P> &M, char transa, char transb) const {
    return gemm (*this, M, transa, transb);
}
template<class T, paradigm P, const bool& b> inline Matrix<T,P>
Matrix<T,P,b>::prodt (const Matrix<T,P> &M) const {
    return gemm (*this, M, 'C');
}
template<class T, paradigm P, const bool& b> inline Matrix<T,P>
Matrix<T,P,b>::operator->* (const Matrix<T,P> &M) const {
    return gemm (*this, M);
}

/**
 * @brief             Matrix vector product A*x
 *
 * Usage: 
 * @code{.cpp}
 *   Matrix<cxdb> A  = rand<cxdb> (20,5);
 *   Matrix<cxdb> x  = rand<cxdb> (20,1);
 *   double prod     = gemv (A, x, 'C');
 * @endcode
 *
 * @param  A          Matrix A
 * @param  x          Vector x
 * @param  trans      Transpose A?
 * 
 * @return            A*x
 */
template<class T> inline Matrix<T> 
gemv (const Matrix<T>& A, const Matrix<T>& x, char trans = 'N') {
    
    assert (isvec(x));
    assert (isvec(A)||is2d(A));
    
	int aw, ah, xh, m, n, one = 1;
	T   alpha = (T)1., beta = (T)0.;
    
	// Column vector
	assert (size(x, 1) == 1);
	
	aw  = (int) size (A, 1);
	ah  = (int) size (A, 0);
	xh  = (int) size (x, 0);
	
	m   = ah;
	n   = aw; 
	
	if (trans == 'N')
		assert (aw == xh);
	else if (trans == 'T' || trans == 'C')
		assert (ah == xh);

	Matrix<T> y ((trans == 'N') ? m : n, 1);
    
	LapackTraits<T>::gemv (trans, m, n, alpha, A.Ptr(), ah, x.Ptr(), one, beta, &y[0], one);
	
	return y;
	
}

/**
 * @brief              Frobenius norm
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl> m = rand<cxfl> (20,10);
 *   float nm       = norm (m); // Lapack driver below produces complex variable with real value
 * @endcode
 *
 * @param  M           Input
 * @return             Eclidean norm
 */
template<class T> inline double
norm (const Matrix<T>& M) {
    
    int n    = (int) numel (M);
    int incx = 1;
    
    return TypeTraits<T>::Real(LapackTraits<T>::nrm2 (n, M.Ptr(), incx));
    
}


/**
 * @brief              Complex dot product (A'*B) on data vector
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<cxdb> a = rand<cxdb> (20,1);
 *   Matrix<cxdb> b = rand<cxdb> (20,1);
 *   cxdb   dotpr   = dotc (a, b);
 * @endcode
 *
 * @param  A           Left factor (is conjugated)
 * @param  B           Right factor
 * @return             A'*B
 */
template <class T> inline T 
dotc (const Matrix<T>& A, const Matrix<T>& B) {
    
	int n = (int) numel(A), one = 1;
	T   res = (T)0.;
    
	assert (n == (int) numel(B));

	LapackTraits<T>::dotc (n, A.Ptr(), one, B.Ptr(), one, &res);
	
	return res;

}


template <class T> inline T 
DOTC (const Matrix<T>& A, const Matrix<T>& B) {
    return dotc (A,B);
}


/**
 * @brief              Dot product (A*B) on data vector
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<float> a = rand<float> (20,1);
 *   Matrix<float> b = rand<float> (20,1);
 *   float dotpr     = dot (a, b);
 * @endcode
 *
 * @param  A           Left factor
 * @param  B           Right factor
 * @return             A*B
 */
template <class T> inline T 
dot  (const Matrix<T>& A, const Matrix<T>& B) {
    
    int n, one;
    T   res;
    
    n   = (int) numel(A);
    assert (n == (int) numel(B));
    
    res = T(0.0);
    one = 1;
    
    LapackTraits<T>::dot (n, A.Ptr(), one, B.Ptr(), one, &res);
    
    return res;
    
}



template <class T> inline T 
DOT  (const Matrix<T>& A, const Matrix<T>& B) {
    return dot (A, B);
}
template<class T, paradigm P, const bool& b> inline T
Matrix<T,P,b>::dotc (const Matrix<T,P>& M) const  {
    return DOTC (*this, M);
}
template<class T, paradigm P, const bool& b> inline T
Matrix<T,P,b>::dot (const Matrix<T,P>& M) const {
    return DOT  (*this, M);
}



#endif // __LAPACK_HPP__
