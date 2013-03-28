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

#include "LapackTraits.hpp"
#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"


/**
 * @brief         Eigenvalue decomposition
 * 
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl> m = rand<cxfl> (10,10), lv, rv;
 *   Matrix<float> ev;
 * 
 *   int res = eig (m, ev, lv, rv, 'N', 'N');
 * @endcode
 * where m is the complex decomposed matrix, lv and rv are the left hand and right 
 * hand side eigen vectors and ev is the real eigenvalue vector.
 *
 * @see           LAPACK driver xGEEV
 *
 * @param  m      Matrix for decomposition
 * @param  ev     Eigenvalues
 * @param  lv     Left  Eigen-vectors
 * @param  rv     Right Eigen-vectors
 * @param  jobvl  Compute left vectors ('N'/'V')
 * @param  jobvr  Compute right vectors ('N'/'V')
 * @return        Status of driver
 */
template <class T, class S> static int
eig (const Matrix<T>& m, Matrix<S>& ev, Matrix<T>& lv, Matrix<T>& rv, const char& jobvl = 'N', const char& jobvr = 'N') {

	if (jobvl != 'N' && jobvl !='V') {
		printf ("EIG Error: Parameter jobvl ('%c' provided) must be 'N' or 'V' \n", jobvl);
		return -1;
	}
	
	if (jobvr != 'N' && jobvr !='V') {
		printf ("EIG Error: Parameter jobvl ('%c' provided) must be 'N' or 'V' \n", jobvr);
		return -1;
	}
	
	// 2D 
	if (!Is2D(m)) {
		printf ("EIG Error: Parameter m must be 2D");
		return -2;
	}
	
	// Square matrix
	if (m.Width() != m.Height()){
		printf ("EIG Error: Parameter m must be square");
		return -3;
	}
	
	int    N     =  size(m, 0);
	
	int    lda   =  N;
	int    ldvl  =  (jobvl) ? N : 1;
	int    ldvr  =  (jobvr) ? N : 1;
	int    info  =  0;
	int    lwork = -1;

    // Appropritely resize the output
    resize (ev,N,1);
    if (jobvl == 'V')
        resize (lv,N,N);
    if (jobvr == 'V')
        resize(rv,N,N);
    

    T* rwork = (T*) malloc (sizeof(T));
	if (typeid(T) == typeid(cxfl) || typeid(T) == typeid(cxdb)) { // Real workspace for complex matrices
		free (rwork);
		rwork = (T*) malloc (N * sizeof(T));
	}
	
	T  wkopt;
	
	Matrix<T> a = m;

	// Work space query
	LapackTraits<T>::geev (&jobvl, &jobvr, N, &a[0], &lda, &ev[0], &lv[0], &ldvl, &rv[0], &ldvr, &wkopt, &lwork, rwork, &info);

	// Initialize work space
	lwork    = (int) creal (wkopt);
	T* work  = (T*) malloc (lwork * sizeof(T));

	// Actual Eigenvalue computation
	LapackTraits<T>::geev (&jobvl, &jobvr, N, &a[0], &lda, &ev[0], &lv[0], &ldvl, &rv[0], &ldvr, &wkopt, &lwork, rwork, &info);

	// Clean up
	free (rwork);
	free (work);
	
	if (info > 0)
		printf ("\nERROR - XGEEV: the QR algorithm failed to compute all the\n eigenvalues, and no eigenvectors have been computed;\n elements %d+1:N of ev contain eigenvalues which\n have converged.\n\n", info) ;
	else if (info < 0)
		printf ("\nERROR - XGEEV: the %d-th argument had an illegal value.\n\n", -info);

	return info;
	
}


/**
 * @brief           Singular value decomposition
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl>  m = rand<cxfl> (20,10), u, v;
 *   Matrix<float> s;
 * 
 *   int res = svd (m, s, u, v, 'N');
 * @endcode
 * where m is the complex decomposed matrix, u and v are the left hand and right 
 * hand side singular vectors and s is the vector of real singular values.
 *
 * @see             LAPACK driver xGESDD
 * 
 * @param  IN       Incoming matrix
 * @param  s        Sorted singular values 
 * @param  U        Left-side singlar vectors
 * @param  V        Right-side singular vectors
 * @param  jobz     Computation mode<br/>
 *                  'A': all m columns of U and all n rows of VT are returned in the arrays u and vt<br/>
 *                  'S', the first min(m, n) columns of U and the first min(m, n) rows of VT are returned in the arrays u and vt;<br/>
 *                  'O', then<br/>
 *                  &nbsp;&nbsp;&nbsp;&nbsp;if m >= n, the first n columns of U are overwritten in the array a and all rows of VT are returned in the array vt;<br/>
 *                  &nbsp;&nbsp;&nbsp;&nbsp;if m < n, all columns of U are returned in the array u and the first m rows of VT are overwritten in the array a;<br/>
 *                  'N', no columns of U or rows of VT are computed.
 * @return          Status of the driver
 */

template<class T, class S> static int 
svd (const Matrix<T>& IN, Matrix<S>& s, Matrix<T>& U, Matrix<T>& V, const char& jobz = 'N') {

    assert (Is2D(IN));
    
	Matrix<T> A (IN);
	
	// SVD only defined on 2D data
	if (!Is2D(A))
		return -2;
	
	int   m, n, lwork, info, lda, mn, ldu = 1, ucol = 1, ldvt = 1, vcol = 1;
	T     wopt;
	void* rwork = 0;
	
	m     =  A.Height();
	n     =  A.Width();
	lwork = -1;
	info  =  0;
	lda   =  m;
	mn    =  MIN(m,n);
	
	if    (((jobz == 'A' || jobz == 'O') && m < n)) 
		ucol = m;
	else if (jobz == 'S')
		ucol = mn;
	
	if      (jobz == 'S' || jobz == 'A' || (jobz == 'O' &&  m < n)) 
		ldu = m;
	
	if      (jobz == 'A' || (jobz == 'O' && m >= n))
		ldvt = n;
	else if (jobz == 'S')
		ldvt = mn;
	
	if      (jobz != 'N')
		vcol = n;
	
	s = resize (s,   mn,    1);
	U = resize (U,  ldu, ucol);
	V = resize (V, ldvt, vcol);
	
	int*   iwork =   (int*) malloc (8 * mn * sizeof(int));
	
	// Only needed for complex data
	if (typeid(T) == typeid(cxfl) || typeid(T) == typeid(cxdb)) {
		if (jobz == 'N') rwork = malloc (mn * 7            * sizeof(T) / 2);
		else             rwork = malloc (mn * (5 * mn + 7) * sizeof(T) / 2);
	}
	
	// Workspace query
	LapackTraits<T>::gesdd (&jobz, &m, &n, &A[0], &lda, &s[0], &U[0], &ldu, &V[0], &ldvt, &wopt, &lwork, rwork, iwork, &info);
	
	// Resize work according to ws query
	lwork   = (int) creal (wopt);
	T* work = (T*) malloc (lwork * sizeof(T));
	
	//SVD
	LapackTraits<T>::gesdd (&jobz, &m, &n, &A[0], &lda, &s[0], &U[0], &ldu, &V[0], &ldvt, work, &lwork, rwork, iwork, &info);

	free (work);
	free (iwork);
	
    //Traspose the baby
	V = !V;
	
	// Clean up
	if (typeid (T) == typeid (cxfl) || typeid (T) == typeid (cxdb)) 
		free (rwork);
	
	if (info > 0)
		printf ("\nERROR - XGESDD: The updating process of SBDSDC did not converge.\n\n");
	else if (info < 0)
		printf ("\nERROR - XGESDD: The %i-th argument had an illegal value.\n\n", -info); 
	
	return info;
	
}
// Convenience calls (for s = svd (A))	
static Matrix<float>
svd (const Matrix<cxfl>& A) {
    Matrix<float> s;
    Matrix<cxfl> u,v;
    svd(A, s, u, v);
    return s;
} 
static Matrix<float>
svd (const Matrix<float>& A) {
    Matrix<float> s;
    Matrix<float> u,v;
    svd(A, s, u, v);
    return s;
} 
static Matrix<double>
svd (const Matrix<cxdb>& A) {
    Matrix<double> s;
    Matrix<cxdb> u,v;
    svd(A, s, u, v);
    return s;
} 
static Matrix<double>
svd (const Matrix<double>& A) {
    Matrix<double> s;
    Matrix<double> u,v;
    svd(A, s, u, v);
    return s;
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
template <class T> static Matrix<T> 
inv (const Matrix<T>& m) {
	
	// 2D 
    assert(Is2D(m));
	
	// Square matrix
	if (size(m,1) != size (m,0)) printf ("Inv Error: Parameter m must be square");
		
	int N = (int) size (m,0);	
	Matrix<T> res = m;
	int  info = 0;
	int *ipiv = (int*) malloc (N * sizeof(int));
	
	// LU Factorisation -------------------

	LapackTraits<T>::getrf (&N, &N, &res[0], &N, ipiv, &info);
	// ------------------------------------
	
	if (info < 0)
		printf ("\nERROR - DPOTRI: the %i-th argument had an illegal value.\n\n", -info);
	else if (info > 1)
		printf ("\nERROR - DPOTRI: the (%i,%i) element of the factor U or L is\n zero, and the inverse could not be computed.\n\n", info, info);
	
	int lwork = -1; 
	T   wopt;
	
	// Workspace determination ------------
	
	LapackTraits<T>::getri (&N, &res[0], &N, ipiv, &wopt, &lwork, &info);
	// ------------------------------------
	
	// Work memory allocation -------------
	
	lwork   = (int) creal (wopt);
	T* work = (T*) malloc (lwork * sizeof(T));
	// ------------------------------------
	
	// Inversion --------------------------
	
	LapackTraits<T>::getri (&N, &res[0], &N, ipiv, work, &lwork, &info);
	// ------------------------------------
	
	free (ipiv);
	free (work);
	
	if (info < 0)
		printf ("\nERROR - XGETRI: The %i-th argument had an illegal value.\n\n", -info);
	else if (info > 0)
		printf ("\nERROR - XGETRI: The leading minor of order %i is not\n positive definite, and the factorization could not be\n completed.", info);
	
	return res;
	
} 


/**
 * @brief                Pseudo invert though SVD
 * 
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl> m = rand<cxfl> (20,10);
 *
 *   m = pinv (m);
 * @endcode
 *
 * @see                  LAPACK driver xGELSD
 * 
 * @param  m             Matrix
 * @param  rcond         Condition number
 * @return               Moore-Penrose pseudoinverse
*/
template<class T> static Matrix<T> 
pinv (const Matrix<T>& m, double rcond = 1.0) {

    assert (Is2D(m));
    
  void *s = 0, *rwork = 0;
	T    *work = 0, wopt = T(0), rwopt = T(0);
	int  *iwork = 0, iwopt = 0;
	
	int  M      =  size(m, 0);
	int  N      =  size(m, 1);
	
	int  nrhs   =  M;
	int  lda    =  M;
	int  ldb    =  MAX(M,N);
	int  lwork  = -1; 
	int  rank   =  0;
	int  info   =  0;
	int  swork  =  sizeof(T) * MIN(M,N);
	
	if (typeid (T) == typeid(cxfl) || typeid (T) == typeid(cxdb))
		swork /= 2;
	
	s      =        malloc (swork);
	
	Matrix<T> b = eye<T>(ldb);

	LapackTraits<T>::gelsd (&M, &N, &nrhs, m.Data(), &lda, &b[0], &ldb, s, rcond, &rank, &wopt, &lwork, &rwopt, &iwopt, &info);
	
	lwork = (int) creal(wopt);
	
	if      (typeid(T) == typeid(cxfl) || 
			 typeid(T) == typeid(cxdb))
		rwork =    malloc ((sizeof(T)/2) * (int) creal(rwopt));
	
	iwork = (int*) malloc (sizeof(int) * iwopt);
	work  = (T*)   malloc (sizeof(T)   * lwork);
	
	LapackTraits<T>::gelsd (&M, &N, &nrhs, m.Data(), &lda, &b[0], &ldb, s, rcond, &rank, work, &lwork,   rwork,  iwork, &info);
	
	if (M > N)
		for (int i = 0; i < M; i++)
			memcpy (&b[i*N], &b[i*M], N * sizeof(T));
	
	b = resize (b, N, M);
	
	if      (typeid (T) == typeid(cxfl) || typeid (T) == typeid(cxdb))
		free (rwork);
	
	free (s);
	free (work);
	free (iwork);
	
	if (info > 0)
		printf ("ERROR XGELSD: the algorithm for computing the SVD failed to converge;\n %i off-diagonal elements of an intermediate bidiagonal form\n did not converge to zero.", info);
	else if (info < 0)
		printf ("ERROR XGELSD: the %i-th argument had an illegal value.", -info);
	
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
template<class T> static Matrix<T> 
chol (const Matrix<T>& A, const char uplo = 'U') {

    assert(Is2D(A));
	
	Matrix<T> res  = A;
	int       info = 0, n = A.Height();
	
	LapackTraits<T>::potrf (&uplo, &n, &res[0], &n, &info);
	
	if (info > 0)
		printf ("\nERROR - XPOTRF: the leading minor of order %i is not\n positive definite, and the factorization could not be\n completed!\n\n", info);
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
template<class T> static Matrix<T> 
gemm (const Matrix<T>& A, const Matrix<T>& B, char transa = 'N', char transb = 'N') {

    assert (Is1D(A)||Is2D(A));
    assert (Is1D(B)||Is2D(B));
    
	int aw, ah, bw, bh, m, n, k, ldc;
	T   alpha, beta;

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
	
	ldc = m;
	
	alpha =  T(1.0);
	beta  =  T(0.0);
	
	Matrix<T> C(m,n);
	
	LapackTraits<T>::gemm (&transa, &transb, &m, &n, &k, &alpha, A.Data(), &ah, B.Data(), &bh, &beta, &C[0], &ldc);

	return C;
	
}



/**
 * @brief              Frobenius norm
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<cxfl> m = rand<cxfl> (20,10);
 *   float normm    = creal(norm (m)); // Lapack driver below produces complex variable with real value
 * @endcode
 *
 * @param  M           Input
 * @return             Eclidean norm
 */
template<class T> static double
norm (const Matrix<T>& M) {

	int n    = (int) numel (M);
	int incx = 1;

	return creal(LapackTraits<T>::nrm2 (n, M.Data(), incx));

}


/**
 * @brief              Complex dot product (A'*B) on data vector
 *
 * Usage:
 * @code{.cpp}
 *   Matrix<cxdb> a = rand<cxdb> (20,1);
 *   Matrix<cxdb> b = rand<cxdb> (20,1);
 *   double dotpr   = dotc (a, b);
 * @endcode
 *
 * @param  A           Left factor (is conjugated)
 * @param  B           Right factor
 * @return             A'*B
 */
template <class T> static T 
dotc (const Matrix<T>& A, const Matrix<T>& B) {

	int n, one;
	T   res;

	n   = (int) numel(A);
	assert (n == (int) numel(B));
	
	res = T(0.0);
	one = 1;
	
	LapackTraits<T>::dotc (n, A.Data(), one, B.Data(), one, &res);
	
	return res;
	
}


template <class T> static T 
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
 *   double dotpr    = dot (a, b);
 * @endcode
 *
 * @param  A           Left factor
 * @param  B           Right factor
 * @return             A*B
 */
template <class T> T 
dot  (const Matrix<T>& A, const Matrix<T>& B) {

	int n, one;
	T   res;

	n   = (int) numel(A);
	assert (n == (int) numel(B));
	
	res = T(0.0);
	one = 1;

	LapackTraits<T>::dot (n, A.Data(), one, B.Data(), one, &res);
	
	return res;
	
}



template <class T> T 
DOT  (const Matrix<T>& A, const Matrix<T>& B) {
	return dot (A, B);
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
 * @param  A          left factor matrix
 * @param  x          Right factor vector
 * @param  trans      Transpose A?
 * 
 * @return            A*x
 */
template<class T> Matrix<T> 
gemv (const Matrix<T>& A, const Matrix<T>& x, char trans = 'N') {

    assert (Is1D(x));
    assert (Is1D(A)||Is2D(A));
    
	int aw, ah, xh, m, n, one;
	T   alpha, beta;

	// Column vector
	assert (size(x, 1) == 1);
	
	aw  = (int) size (A, 1);
	ah  = (int) size (A, 0);
	xh  = (int) size (x, 0);
	
	m   = ah;
	n   = aw; 
	one = 1;
	
	if (trans == 'N')
		assert (aw == xh);
	else if (trans == 'T' || trans == 'C')
		assert (ah == xh);
	
	alpha  = T(1.0);
	beta   = T(0.0);
	
	Matrix<T> y ((trans == 'N') ? m : n, 1);

	LapackTraits<T>::gemv (&trans, &m, &n, &alpha, A.Data(), &ah, x.Data(), &one, &beta, &y[0], &one);
	
	return y;
	
}



#endif // __LAPACK_HPP__
