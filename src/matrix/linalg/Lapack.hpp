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
#include "Symmetry.hpp"
#include "CX.hpp"

#ifdef HAVE_CXX11_TUPLE
#include <tuple>
#define TUPLE std::tuple
#define GET std::get
#else
#include <boost/tuple/tuple.hpp>
#define TUPLE boost::tuple
#define GET boost::get
#endif


template<class T> struct eig_t {
	Matrix<T> lv;
	Matrix<typename TypeTraits<T>::CT> ev;
	Matrix<T> rv;
};

template<class T> inline static eig_t<T> eigs (const Matrix<T>& A, char jobz = 'V') {
	typedef typename TypeTraits<T>::RT real;
	typedef typename TypeTraits<T>::CT cplx;
	char uplo = 'U';
	assert((TypeTraits<T>::IsReal()&&issymmetric(A)) || (TypeTraits<T>::IsComplex()&&ishermitian(A)));
	int n = size(A,0), lwork = -1, liwork = -1, lrwork=-1, info;
	eig_t<T> e;
	Vector<real> rwork(1), w(n);
	Vector<int> iwork(1);
	Vector<T> work(1);
	e.lv = A;
	e.ev = Matrix<cplx>(n,1);
	LapackTraits<T>::syevd (jobz, uplo, n, e.lv.Container(), w, work, lwork, rwork, lrwork, iwork, liwork, info);
	lwork = 1.5*TypeTraits<T>::Real(work[0]); work.resize(lwork);
	liwork = iwork[0]; iwork.resize(liwork);
	lrwork = (TypeTraits<T>::IsComplex()) ? rwork[0] : 1;
	rwork.resize(lrwork);
	LapackTraits<T>::syevd (jobz, uplo, n, e.lv.Container(), w, work, lwork, rwork, lrwork, iwork, liwork, info);
    if (info > 0) {
        printf ("\nERROR**: X(SY/HE)VD: ");
        if (jobz=='N') {
        	printf ("the algorithm failed to converge; %dth off-diagonal elements of an intermediate "
        			"tridiagonal form did not converge to zero;\n\n", info);
        } else {
        	printf ("the algorithm failed to compute an eigenvalue while working on the submatrix "
        			"lying in rows and columns %d through %d.", info/(n+1), info%(n+1));

        }
    } else if (info < 0) {
        printf ("\nERROR - XGEEV: the %d-th argument had an illegal value.\n\n", -info);
    }
	for (size_t i = 0; i < n; ++i)
		e.ev[i] = w[i];
	return e;
}


/**
 * @brief         Eigenvalue decomposition
 * 
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl> m = rand<cxfl> (10,10), ev;
 *   Matrix<float> lv, rv;
 *   TUPLE<Matrix<float>,Matrix<cxfl>,Matrix<float>> ler
 *   ler = eig2 (m);
 *   lv = GET<0>(ler)
 *   ev = GET<1>(ler);
 *   rv = GET<2>(ler)
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
template <class T> inline static eig_t<T> eigu (const Matrix<T>& m, char jobvl = 'V', char jobvr = 'N') {
    
    typedef typename TypeTraits<T>::CT CT;
    typedef typename TypeTraits<T>::RT RT;
    eig_t<T> ret;
    
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
    ret.ev = Matrix<CT>(N,1);
    
    if (jobvl == 'V') ret.lv = Matrix<T>(N,N);
    if (jobvr == 'V') ret.rv = Matrix<T>(N,N);
    
    // Workspace for real numbers
    Vector<RT> rwork = (TypeTraits<T>::IsComplex()) ?
    		Vector<RT>(2*N) : Vector<RT>(1);

    // Workspace for complex numbers
    Vector<T>  work = Vector<T>(1);
    
    // Need copy. Lapack destroys A on output.
    Matrix<T> a = m;

    // Work space query
    LapackTraits<T>::geev (jobvl, jobvr, N, a.Ptr(), lda, ret.ev.Ptr(), ret.lv.Ptr(),
    		ldvl, ret.rv.Ptr(), ldvr, work.ptr(), lwork, rwork.ptr(), info);
    
    // Initialize work space
    lwork = (int) TypeTraits<T>::Real (work[0]);
    work.resize(lwork);
    
    // Actual Eigenvalue computation
    LapackTraits<T>::geev (jobvl, jobvr, N, a.Ptr(), lda, ret.ev.Ptr(), ret.lv.Ptr(),
    		ldvl, ret.rv.Ptr(), ldvr, work.ptr(), lwork, rwork.ptr(), info);

    if (info > 0)
        printf ("\nERROR - XGEEV: the QR algorithm failed to compute all the\n "
        		"eigenvalues, and no eigenvectors have " \
                "been computed;\n elements %d+1:N of ev contain eigenvalues "
                "which\n have converged.\n\n", info) ;
    else if (info < 0)
        printf ("\nERROR - XGEEV: the %d-th argument had an illegal value.\n\n", -info);
    
    return ret;
    
}
template<class T> inline static eig_t<T> eig2 (const Matrix<T>& A, char c1 = 'V', char c2 = 'N') {
	if ((TypeTraits<T>::IsReal()&&issymmetric(A)) || (TypeTraits<T>::IsComplex()&&ishermitian(A)))
		return eigs (A, c1);
	else
		return eigu (A, c1, c2);
}
template<class T> inline static Matrix<typename TypeTraits<T>::CT> eig (const Matrix<T>& m) {
	eig_t<T> e = eig2<T>(m, 'N', 'N');
	return e.ev;
}



/**
 * @brief           Singular value decomposition
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl>  m = rand<cxfl> (20,10), u, v;
 *   Matrix<float> s;
 *   TUPLE<Matrix<cxfl>, Matrix<float>, Matrix<cxfl>> usv;
 *   usv = svd2 (m);
 *   u = GET<0>(usv);
 *   s = GET<1>(usv);
 *   v = GET<2>(usv);
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

template<class T> inline TUPLE<Matrix<T>,Matrix<typename TypeTraits<T>::RT>,Matrix<T> >
svd2 (const Matrix<T>& M, char jobz = 'N') {
    
    typedef typename TypeTraits<T>::RT RT;
    TUPLE<Matrix<T>,Matrix<RT>,Matrix<T> > ret;
    
    assert (is2d(M));
    assert (jobz == 'N' || jobz == 'S' || jobz == 'A');
    
    Matrix<T> A (M);
    
    int   m, n, lwork, info = 0, lda, mn, ldu = 1, ucol = 1, ldvt = 1, vtcol = 1;

    Vector<RT> rwork;
    Vector<T>  work = Vector<T>(1);
    
    m     =  A.Height(); n = A.Width();
    lwork = -1;
    lda   =  m;
    mn    =  std::min(m,n);
    
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

    Matrix<RT>& s = GET<1>(ret) = Matrix<RT>(  mn,  IONE);
    Matrix<T>&  U = GET<0>(ret) = Matrix<T> ( ldu,  ucol);
    Matrix<T>&  V = GET<2>(ret) = Matrix<T> (ldvt, vtcol);
    
    Vector<int> iwork = Vector<int>(8 * mn);
    
    size_t nr = (typeid(T) == typeid(cxfl) || typeid(T) == typeid(cxdb)) ?
            ((jobz == 'N') ? mn * 7 : mn * (5 * mn + 7)) : 1;
    rwork.resize(nr);

    // Workspace query
    LapackTraits<T>::gesdd (jobz, m, n, A.Ptr(), lda, s.Ptr(), U.Ptr(), ldu, V.Ptr(),
            ldvt, work.ptr(), lwork, rwork.ptr(), iwork.ptr(), info);
    
    lwork = (int) TypeTraits<T>::Real(work[0]);
    work.resize(lwork);
    
    // SVD
    LapackTraits<T>::gesdd (jobz, m, n, A.Ptr(), lda, s.Ptr(), U.Ptr(), ldu, V.Ptr(),
            ldvt, work.ptr(), lwork, rwork.ptr(), iwork.ptr(), info);
    
    
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
template<class T>
inline Matrix<typename TypeTraits<T>::RT> svd (const Matrix<T>& A) {
    return GET<1>(svd2<T>(A));
}
template<class T> inline Matrix<T> pca (const Matrix<T>& A) {
	return GET<2>(svd2<T>(A));
}
template<class T>
inline TUPLE<Matrix<typename TypeTraits<T>::RT>,Matrix<T>> pca2 (const Matrix<T>& A) {
	typedef typename TypeTraits<T>::RT RT;
	TUPLE<Matrix<T>,Matrix<typename TypeTraits<T>::RT>,Matrix<T> > ret = svd2<T>(A,'S');
	Matrix<RT>& S = GET<1>(ret);
	S *= S/numel(S);
	Matrix<T>& V = GET<2>(ret);
	return TUPLE<Matrix<typename TypeTraits<T>::RT>,Matrix<T> > (S,V);
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
template <class T> inline Matrix<T> inv (const Matrix<T>& m) {
    
    // 2D 
    assert(issquare(m));
    
    int N = (int) size (m,0), lwork = -1, info = 0;
    Matrix<T> res = m;
    Vector<int> ipiv = Vector<int>(N);
    
    // LU Factorisation -------------------
    LapackTraits<T>::getrf (N, N, res.Ptr(), N, ipiv.ptr(), info);
    
    if (info < 0)
        printf ("\nERROR - DPOTRI: the %i-th argument had an illegal value.\n\n", -info);
    else if (info > 1)
        printf ("\nERROR - DPOTRI: the (%i,%i) element of the factor U or L is\n "
        		"zero, and the inverse could not be computed.\n\n", info, info);
    
    // Workspace determination ------------
    Vector<T> work = Vector<T>(1);
    LapackTraits<T>::getri (N, res.Ptr(), N, ipiv.ptr(), work.ptr(), lwork, info);
    
    // Inversion --------------------------
    lwork = (int) TypeTraits<T>::Real (work[0]);
    work.resize(lwork);
    LapackTraits<T>::getri (N, res.Ptr(), N, ipiv.ptr(), work.ptr(), lwork, info);

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
template<class T> inline Matrix<T> pinv (const Matrix<T>& m, char trans = 'N') {
    
    assert (is2d(m));
    
    int  M      =  size(m, 0);
    int  N      =  size(m, 1);
    int  nrhs   =  M;
    int  lda    =  M;
    int  ldb    =  std::max(M,N);
    int  lwork  = -1;
    int  info   =  0;

    Matrix<T> mm (m), b = eye<T>(ldb);
    Vector<T> work = Vector<T>(1);
    
    // Workspace determination
    LapackTraits<T>::gels (trans, M, N, nrhs, mm, lda, b, ldb, work, lwork, info);
    
    // Solving
    lwork = (int) TypeTraits<T>::Real(work[0]);
    work.resize(lwork);
    LapackTraits<T>::gels (trans, M, N, nrhs, mm, lda, b, ldb, work, lwork, info);
    
    size_t cpsz = N * sizeof(T);
    if (M > N)
        for (int i = 0; i < M; i++)
            memcpy (&b[i*N], &b[i*M], cpsz);
    
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
    
    LapackTraits<T>::potrf (uplo, n, res.Ptr(), n, info);
    
    if (info > 0)
        printf ("\nERROR - XPOTRF: the leading minor of order %i is not\n positive definite, and the factorization " \
                "could not be\n completed!\n\n", info);
    else if (info < 0)
        printf ("\nERROR - XPOTRF: the %i-th argument had an illegal value.\n\n!", -info);
    
    for (size_t i = 0; i < (size_t)n-1; i++)
        for (size_t j = i+1; j < (size_t)n; j++)
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
	
	LapackTraits<T>::gemm (transa, transb, m, n, k, alpha, A.Ptr(), ah, B.Ptr(), bh, beta, C.Ptr(), m);
    
	return C;
	
}

template<class T, paradigm P> Matrix<T,P>
Matrix<T,P>::prod  (const Matrix<T,P> &M, char transa, char transb) const {
    return gemm (*this, M, transa, transb);
}
template<class T, paradigm P> inline Matrix<T,P>
Matrix<T,P>::prodt (const Matrix<T,P> &M) const {
    return gemm (*this, M, 'C');
}
template<class T, paradigm P> inline Matrix<T,P>
Matrix<T,P>::operator->* (const Matrix<T,P> &M) const {
    return gemm (*this, M);
}
template<class T> Matrix<std::complex<T> >
gemm (const Matrix<T>& A, const Matrix<std::complex<T> >& B, char transa = 'N', char transb = 'N') {
	Matrix<std::complex<T> > AC = A;
	return gemm(AC,B,transa,transb);
}
template<class T> Matrix<std::complex<T> >
gemm (const Matrix<std::complex<T> >& A, const Matrix<T>& B, char transa = 'N', char transb = 'N') {
	Matrix<std::complex<T> > BC = B;
	return gemm(A,BC,transa,transb);
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
    
	LapackTraits<T>::gemv (trans, m, n, alpha, A.Ptr(), ah, x.Ptr(), one, beta, y.Ptr(), one);
	
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
 * @param  what        Which norm (E, F, 0, 1)
 * @return             Eclidean norm ('E', default) 0, Frobenius, 1 ('0', 'F', '1')
 */
template<class T> inline typename LapackTraits<T>::RType
norm (const Matrix<T>& M, const char what = 'E') {
	switch (what)
	{
		case 'E': return LapackTraits<T>::nrm2 (M.Size(), M.Ptr(), 1); break;
		case 'e': return LapackTraits<T>::nrm2 (M.Size(), M.Ptr(), 1); break;
		default : return LapackTraits<T>::lange (what, size(M,0), size(M,1), M.Container()); break;
	}
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
	T res;// = (T)0.;
    
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
    
    int n, one = 1;
    T   res;
    
    n   = (int) numel(A);
    assert (n == (int) numel(B));
    
    LapackTraits<T>::dot (n, A.Ptr(), one, B.Ptr(), one, &res);
    
    return res;
    
}



template <class T> inline T 
DOT  (const Matrix<T>& A, const Matrix<T>& B) {
    return dot (A, B);
}
template<class T, paradigm P> inline T
Matrix<T,P>::dotc (const Matrix<T,P>& M) const  {
    return DOTC (*this, M);
}
template<class T, paradigm P> inline T
Matrix<T,P>::dot (const Matrix<T,P>& M) const {
    return DOT  (*this, M);
}



#endif // __LAPACK_HPP__
