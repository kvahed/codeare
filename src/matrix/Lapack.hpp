
/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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

extern "C" {

	// Cholesky factorization of a complex Hermitian positive definite matrix
	void cpotrf_ (char* uplo, int* n, void* a, int* lda, int *info);
	void dpotrf_ (char* uplo, int* n, void* a, int* lda, int *info);
	
	// Computes an LU factorization of a general M-by-N matrix A
	void cgetrf_ (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	void dgetrf_ (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	void zgetrf_ (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	void sgetrf_ (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	
	// Inverse of a complex Hermitian positive definite matrix using cpotrf/cpptrf
	void cpotri_ (char* uplo, int*n, void *a, int* lda, int*info);
	void dpotri_ (char* uplo, int*n, void *a, int* lda, int*info);
	
	// Matrix inversion through cholesky decomposition
	void cgetri_ (int *n, void *a, int* lda, int *ipiv, void *work, int *lwork, int*info);
	void dgetri_ (int *n, void *a, int* lda, int *ipiv, void *work, int *lwork, int*info);
	void zgetri_ (int *n, void *a, int* lda, int *ipiv, void *work, int *lwork, int*info);
	void sgetri_ (int *n, void *a, int* lda, int *ipiv, void *work, int *lwork, int*info);
	
	// Eigen value computations
	void cgeev_  (char *jvl, char *jvr, int *n, void *a, int *lda, void *w ,           void *vl, int *ldvl, void *vr, int *ldvr, void *work, int *lwork, void *rwork, int *info);
	void zgeev_  (char *jvl, char *jvr, int *n, void *a, int *lda, void *w ,           void *vl, int *ldvl, void *vr, int *ldvr, void *work, int *lwork, void *rwork, int *info);
	void dgeev_  (char *jvl, char *jvr, int *n, void *a, int *lda, void *wr, void *wi, void *vl, int *ldvl, void *vr, int *ldvr, void *work, int *lwork,              int *info);
	void sgeev_  (char *jvl, char *jvr, int *n, void *a, int *lda, void *wr, void *wi, void *vl, int *ldvl, void *vr, int *ldvr, void *work, int *lwork,              int *info);
	
	// Singular value decomposition 
	void cgesdd_ (const char *jobz, int*m, int *n, void *a, int *lda, void *s, void*u, int*ldu, void *vt, int *ldvt, void *work, int*lwork, void *rwork, int *iwork, int*info);
	void zgesdd_ (const char *jobz, int*m, int *n, void *a, int *lda, void *s, void*u, int*ldu, void *vt, int *ldvt, void *work, int*lwork, void *rwork, int *iwork, int*info);
	void dgesdd_ (const char *jobz, int*m, int *n, void *a, int *lda, void *s, void*u, int*ldu, void *vt, int *ldvt, void *work, int*lwork,               int *iwork, int*info);
	void sgesdd_ (const char *jobz, int*m, int *n, void *a, int *lda, void *s, void*u, int*ldu, void *vt, int *ldvt, void *work, int*lwork,               int *iwork, int*info);
	
	// Pseudo-inversion 
	void zgelsd_ (int* m, int* n, int* nrhs, void* a, int* lda, void* b, int* ldb, void* s, void* rcond, int* rank, void* work, int* lwork, void* rwork, int* iwork, int* info);
	void cgelsd_ (int* m, int* n, int* nrhs, void* a, int* lda, void* b, int* ldb, void* s, void* rcond, int* rank, void* work, int* lwork, void* rwork, int* iwork, int* info);
	void dgelsd_ (int* m, int* n, int* nrhs, void* a, int* lda, void* b, int* ldb, void *s, void* rcond, int* rank, void* work, int* lwork, void*             iwork, int* info);
	void sgelsd_ (int* m, int* n, int* nrhs, void* a, int* lda, void* b, int* ldb, void *s, void* rcond, int* rank, void* work, int* lwork, void*             iwork, int* info);

	// Matrix vector multiplication
	void sgemv_  (char* trans, int* m, int* n, void* alpha, void *a, int* lda, void *x, int* incx, void* beta, void *y, int* incy);
	void dgemv_  (char* trans, int* m, int* n, void* alpha, void *a, int* lda, void *x, int* incx, void* beta, void *y, int* incy);
	void cgemv_  (char* trans, int* m, int* n, void* alpha, void *a, int* lda, void *x, int* incx, void* beta, void *y, int* incy);
	void zgemv_  (char* trans, int* m, int* n, void* alpha, void *a, int* lda, void *x, int* incx, void* beta, void *y, int* incy);

	// Matrix matrix multiplication
	void sgemm_  (char *transa, char *transb, int  *m, int   *n, int *k, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
	void dgemm_  (char *transa, char *transb, int  *m, int   *n, int *k, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
	void cgemm_  (char *transa, char *transb, int  *m, int   *n, int *k, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
	void zgemm_  (char *transa, char *transb, int  *m, int   *n, int *k, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
	
}

#include "Matrix.hpp"
#include "cblas.h"

class Lapack {

	public:

	/**
	 * @brief         Eigenvalue decomposition
	 *                Eigenvalue decomposition with Lapack routines Xgeev
	 *
	 * @param  m      Matrix for decomposition
	 * @param  lv     Left  Eigen-vectors
	 * @param  rv     Right Eigen-vectors
	 * @param  jobvl  Compute left vectors ('N'/'V')
	 * @param  jobvr  Compute right vectors ('N'/'V')
	 * @return        Eigen values
	 */
	template <class T, class S> static int
	EIG (Matrix<T>& m, Matrix<S>& ev, Matrix<T>& lv, Matrix<T>& rv, char jobvl = 'N', char jobvr = 'N') {

		if (jobvl != 'N' && jobvl !='V') {
			printf ("EIG Error: Parameter jobvl ('%c' provided) must be 'N' or 'V' \n", jobvl);
			return -1;
		}

		if (jobvr != 'N' && jobvr !='V') {
			printf ("EIG Error: Parameter jobvl ('%c' provided) must be 'N' or 'V' \n", jobvr);
			return -1;
		}

		// 2D 
		if (!m.Is2D()) {
			printf ("EIG Error: Parameter m must be 2D");
			return -2;
		}

		// Square matrix
		if (m.Width() != m.Height()){
			printf ("EIG Error: Parameter m must be square");
			return -3;
		}
		
		int    N     =  m.Height();

		int    lda   =  N;
		int    ldvl  =  (jobvl) ? N : 1;
		int    ldvr  =  (jobvr) ? m.Width() : 1;
		int    info  =  0;
		int    lwork = -1;
		
		T*     w;
		T*     wi;
		T*     rwork;

		if (typeid(T) == typeid(float) || typeid(T) == typeid(double)) {    // Complex eigen values for real matrices
			w     = (T*) malloc (N * sizeof(T));
			wi    = (T*) malloc (N * sizeof(T));
		} else if (typeid(T) == typeid(cxfl) || typeid(T) == typeid(cxdb))  // Real workspace for complex matrices
			rwork = (T*) malloc (N * sizeof(T));
	  
		T  wkopt;

		// Workspace query
		if (typeid(T) == typeid(cxfl))
			cgeev_ (&jobvl, &jobvr, &N, &m[0], &lda, &ev[0], &lv[0], &ldvl, &rv[0], &ldvr, &wkopt, &lwork, rwork, &info);
		else if (typeid(T) == typeid(cxdb))
			zgeev_ (&jobvl, &jobvr, &N, &m[0], &lda, &ev[0], &lv[0], &ldvl, &rv[0], &ldvr, &wkopt, &lwork, rwork, &info);
		else if (typeid(T) == typeid(double))
			dgeev_ (&jobvl, &jobvr, &N, &m[0], &lda,  w, wi, &lv[0], &ldvl, &rv[0], &ldvr, &wkopt, &lwork,        &info);
		else if (typeid(T) == typeid(float))
			sgeev_ (&jobvl, &jobvr, &N, &m[0], &lda,  w, wi, &lv[0], &ldvl, &rv[0], &ldvr, &wkopt, &lwork,        &info);
		
		// Intialise work space
		lwork    = (int) creal (wkopt);
		T* work  = (T*) malloc (lwork * sizeof(T));

		// Actual eigen value comp
		if (typeid(T) == typeid(cxfl)) {
			cgeev_ (&jobvl, &jobvr, &N, &m[0], &lda, &ev[0], &lv[0], &ldvl, &rv[0], &ldvr, work, &lwork, rwork, &info);
		} else if (typeid(T) == typeid(cxdb)) {
			zgeev_ (&jobvl, &jobvr, &N, &m[0], &lda, &ev[0], &lv[0], &ldvl, &rv[0], &ldvr, work, &lwork, rwork, &info);
		} else if (typeid(T) == typeid(double)) {
			dgeev_ (&jobvl, &jobvr, &N, &m[0], &lda,  w, wi, &lv[0], &ldvl, &rv[0], &ldvr, work, &lwork,        &info);
			for (size_t i = 0; i < N; i++) {
				double f[2] = {((double*)w)[i], ((double*)wi)[i]};
				memcpy(&ev[i], f, 2 * sizeof(double));
			}
		} else if (typeid(T) == typeid(float)) {
			sgeev_ (&jobvl, &jobvr, &N, &m[0], &lda,  w, wi, &lv[0], &ldvl, &rv[0], &ldvr, work, &lwork,        &info);
			for (size_t i = 0; i < N; i++) {
				float f[2] = {((float*)w)[i], ((float*)wi)[i]};
				memcpy(&ev[i], f, 2 * sizeof(float));
			}
		}

		
		// Clean up
		if (typeid(T) == typeid(float) || typeid(T) == typeid(double)) {
			free (w);
			free (wi);
		} else if (typeid(T) == typeid(cxfl) || typeid(T) == typeid(cxdb)) 
			free (rwork);
		free (work);

		
		return info;
		
	}
	

	template<class T, class S> static int 
	SVD (Matrix<T>& A, Matrix<S>& s, Matrix<T>& U, Matrix<T>& V, const char jobz = 'N') {
		
		// SVD only defined on 2D data
		if (!A.Is2D())
			return -2;

		int   m, n, lwork, info, lda, mn, ldu = 1, ucol = 1, ldvt = 1, vcol = 1;
		T     wopt;
		void* rwork;

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

		s.Resize (mn,1);
		U.Resize (ldu,ucol);
		V.Resize (ldvt,vcol);

		int*   iwork =   (int*) malloc (8 * mn * sizeof(int));
		
		// Only needed for complex data
		if (typeid(T) == typeid(cxfl) || typeid(T) == typeid(cxdb)) {
			if (jobz == 'N') rwork = malloc (mn * 7            * sizeof(T) / 2);
			else             rwork = malloc (mn * (5 * mn + 7) * sizeof(T) / 2);
		}
		
		// Workspace query
		if      (typeid(T) == typeid(cxfl))
			cgesdd_ (&jobz, &m, &n, &A[0], &lda, &s[0], &U[0], &ldu, &V[0], &ldvt, &wopt, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(cxdb))
			zgesdd_ (&jobz, &m, &n, &A[0], &lda, &s[0], &U[0], &ldu, &V[0], &ldvt, &wopt, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(double))
			dgesdd_ (&jobz, &m, &n, &A[0], &lda, &s[0], &U[0], &ldu, &V[0], &ldvt, &wopt, &lwork,        iwork, &info);
		else if (typeid(T) == typeid(float))
			sgesdd_ (&jobz, &m, &n, &A[0], &lda, &s[0], &U[0], &ldu, &V[0], &ldvt, &wopt, &lwork,        iwork, &info);
		
		// Resize work according to ws query
		lwork   = (int) creal (wopt);
		T* work = (T*) malloc (lwork * sizeof(T));
		
		//SVD
		if      (typeid(T) == typeid(cxfl))
			cgesdd_ (&jobz, &m, &n, &A[0], &lda, &s[0], &U[0], &ldu, &V[0], &ldvt, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(cxdb))
			zgesdd_ (&jobz, &m, &n, &A[0], &lda, &s[0], &U[0], &ldu, &V[0], &ldvt, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(double))
			dgesdd_ (&jobz, &m, &n, &A[0], &lda, &s[0], &U[0], &ldu, &V[0], &ldvt, work, &lwork,        iwork, &info);
		else if (typeid(T) == typeid(float))
			sgesdd_ (&jobz, &m, &n, &A[0], &lda, &s[0], &U[0], &ldu, &V[0], &ldvt, work, &lwork,        iwork, &info);

		// Transpose V
		//V = !V;

		// Clean up
		//if (typeid (T) == typeid (cxfl) || typeid (T) == typeid (cxdb)) free (rwork);

		//free (work);
		//free (iwork);
		
		if (info > 0)
			printf ("The updating process of SBDSDC did not converge.\n");
		else if (info < 0)
			printf ("the %i-th argument had an illegal value.\n", -info); 
		
		return info;
		
	} 
	

	template <class T> static Matrix<T> 
	Inv (Matrix<T>& m) {
		
		// 2D 
		if (!m.Is2D())
			printf ("Inv Error: Parameter m must be 2D");

		// Square matrix
		if (m.Width() != m.Height())
			printf ("Inv Error: Parameter m must be square");
		
		Matrix<T> res = m;
		
		int  N    = m.Height();
		int  info = 0;
		int *ipiv = (int*) malloc (N * sizeof(int));

		char up = 'U';
		char lo = 'L';
		
		// LQ Factorisation -------------------
		
		if (typeid (T) == typeid (cxfl)) 
			cgetrf_ (&N, &N, &res[0], &N, ipiv, &info);
		else if (typeid (T) == typeid (double))
			dgetrf_ (&N, &N, &res[0], &N, ipiv, &info);
		else if (typeid (T) == typeid (cxdb)) 
			zgetrf_ (&N, &N, &res[0], &N, ipiv, &info);
		else if (typeid (T) == typeid (float))
			sgetrf_ (&N, &N, &res[0], &N, ipiv, &info);
		
		// ------------------------------------
		
		int lwork = -1; 
		T*  work  = (T*) malloc (sizeof(T));
		
		// Workspace determination ------------
		
		if (typeid (T) == typeid (cxfl)) 
			cgetri_ (&N, &res[0], &N, ipiv, work, &lwork, &info);
		else if (typeid (T) == typeid (cxdb)) 
			zgetri_ (&N, &res[0], &N, ipiv, work, &lwork, &info);
		else if (typeid (T) == typeid (double))
			dgetri_ (&N, &res[0], &N, ipiv, work, &lwork, &info);
		else if (typeid (T) == typeid (float))
			sgetri_ (&N, &res[0], &N, ipiv, work, &lwork, &info);
		
		// Work memory allocation -------------
	
		if (typeid (T) == typeid(cxfl) || typeid (T) == typeid(float))
			lwork = (int) ((float*)work)[0];
		else if (typeid (T) == typeid(cxdb) || typeid (T) == typeid(double))
			lwork = (int) ((double*)work)[0];

		work  = (T*)  realloc (work, lwork * sizeof(T));

		// Inversion --------------------------

		if (typeid (T) == typeid (cxfl)) 
			cgetri_ (&N, &res[0], &N, ipiv, work, &lwork, &info);
		else if (typeid (T) == typeid (cxdb)) 
			zgetri_ (&N, &res[0], &N, ipiv, work, &lwork, &info);
		else if (typeid (T) == typeid (double))
			dgetri_ (&N, &res[0], &N, ipiv, work, &lwork, &info);
		else if (typeid (T) == typeid (float))
			sgetri_ (&N, &res[0], &N, ipiv, work, &lwork, &info);
		// ------------------------------------
		
		free (ipiv);
		free (work);
		
		return res;
		
	} 
	


	template<class T> static Matrix<T> 
	Pinv (Matrix<T>& m) {
		
		void *s, *rwork;
		T    *work;
		int  *iwork;
		
		int  M      =  m.Height();
		int  N      =  m.Width();

		int  nrhs   =  M;
		int  lda    =  M;
		int  ldb    =  MAX(M,N);
		int  lwork  = -1; 
		int  lrwork =  0;
		int  liwork =  0;
		int  rank   =  0;
		int  info   =  0;
		int  swork  =  sizeof(T) * MIN(M,N);
		
		work   = (T*)   malloc  (sizeof(T)); 
		iwork  = (int*) malloc  (sizeof(int)); 
		rwork  =        malloc  (sizeof(T) / 2);    


		if (typeid (T) == typeid(cxfl) || typeid (T) == typeid(cxdb))
			swork /= 2;

		s      =        malloc (swork);    

		Matrix<T> b =  Matrix<T>::Id(ldb);
		
		double drcond  = -1.0;
		float  frcond  = -1.0;
		
		if (typeid(T) == typeid(cxfl))
			cgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &frcond, &rank, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(cxdb))
			zgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &drcond, &rank, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(double))
			dgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &drcond, &rank, work, &lwork,        iwork, &info);
		else if (typeid(T) == typeid(float))
			sgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &frcond, &rank, work, &lwork,        iwork, &info);
		
		liwork = iwork[0];

		if      (typeid (T) == typeid(cxfl))                                 lrwork = (int) ((float*) rwork)[0];
		else if (typeid (T) == typeid(cxdb))                                 lrwork = (int) ((double*)rwork)[0];

		if      (typeid (T) == typeid(cxfl) || typeid (T) == typeid(float))  lwork  = (int) ((float*) work)[0];
		else if (typeid (T) == typeid(cxdb) || typeid (T) == typeid(double)) lwork  = (int) ((double*)work)[0];


		if      (typeid (T) == typeid(cxfl) || typeid (T) == typeid(cxdb))
			rwork =    malloc ((sizeof(T)/2)*lrwork);
		
		iwork = (int*) malloc (liwork*sizeof(int));
		work  = (T*)   malloc (lwork *sizeof(T));
		
		if (typeid(T) == typeid(cxfl))
			cgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &frcond, &rank, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(cxdb))
			zgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &drcond, &rank, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(double))
			dgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &drcond, &rank, work, &lwork,        iwork, &info);
		else if (typeid(T) == typeid(float))
			sgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &frcond, &rank, work, &lwork,        iwork, &info);

		if (M > N)
			for (size_t i = 0; i < M; i++)
				memcpy (&b[i*N], &b[i*M], N * sizeof(T));

		b.Resize(N,M);
		
		free (s);
		free (work);
		free (rwork);
		free (iwork);
		
		return b;
		
	}
	
	

	template<class T> static Matrix<T> 
	Cholesky (Matrix<T>& m, const char uplo) {
		
		Matrix<T> res = m;
		int       info = 0;
		
		if (typeid(T) == typeid(cxfl))
			//dpotrf_ ( &uplo, &tmp.Width(), &tmp[0], &tmp.Height(), &info);
		
		return res;
		
	}
	

	template<class T> static Matrix<T> 
	GEMM (Matrix<T>& a, Matrix<T>& b, char transa = 'N', char transb = 'N') {
		
		int aw = (int)a.Width(), ah = (int)a.Height(), bw = (int)b.Width(), bh = (int)b.Height();
		
		if      ( transa == 'N'                   &&  transb == 'N'                  ) assert (aw == bh);
		else if ( transa == 'N'                   && (transb == 'T' || transb == 'C')) assert (aw == bw);
		else if ((transa == 'T' || transa == 'C') &&  transb == 'N'                  ) assert (ah == bh);
		else if ((transa == 'T' || transa == 'C') && (transb == 'T' || transb == 'C')) assert (ah == bw);
		
		int m, n, k, ldc, cm, cn;
		
		if (transa == 'N') {
			m   = (int) ah;
			k   = (int) aw;
		} else if (transa == 'T' || transa == 'C') {
			m   = (int) aw;
			k   = (int) ah;
		}

		if (transb == 'N')
			n   = (int) bw;
		else if (transb == 'T' || transb == 'C')
			n   = (int) bh;
		
		ldc = m;

		T    alpha  =       T(1.0);
		T    beta   =       T(0.0);
		
		Matrix<T> c (m, n);
		
		if      (typeid(T) == typeid(double))
			dgemm_ (&transa, &transb, &m, &n, &k, &alpha, &a[0], &ah, &b[0], &bh, &beta, &c[0], &ldc);
		else if (typeid(T) == typeid(float))
			sgemm_ (&transa, &transb, &m, &n, &k, &alpha, &a[0], &ah, &b[0], &bh, &beta, &c[0], &ldc);
		else if (typeid(T) == typeid(cxfl))
			cgemm_ (&transa, &transb, &m, &n, &k, &alpha, &a[0], &ah, &b[0], &bh, &beta, &c[0], &ldc);
		else if (typeid(T) == typeid(cxdb))
			zgemm_ (&transa, &transb, &m, &n, &k, &alpha, &a[0], &ah, &b[0], &bh, &beta, &c[0], &ldc);
		
		return c;
		
	}



	template<class T> static Matrix<T> 
	GEMV (Matrix<T>& A, Matrix<T>& x, char trans = 'N') {

		assert (x.Width() == 1);
		
		int aw = (int) A.Width(), ah = (int) A.Height(), xh = (int) x.Height();

		int m = ah, n = aw, one = 1;
		
		if (trans == 'N')
			assert (aw == xh);
		else if (trans == 'T' || trans == 'C')
			assert (ah == xh);
		
		T    alpha  = T(1.0);
		T    beta   = T(0.0);
		
		Matrix<T> y ((trans == 'N') ? m : n, 1);
		
		if      (typeid(T) == typeid(double)) dgemv_ (&trans, &m, &n, &alpha, &A[0], &ah, &x[0], &one, &beta, &y[0], &one);
		else if (typeid(T) == typeid(float))  sgemv_ (&trans, &m, &n, &alpha, &A[0], &ah, &x[0], &one, &beta, &y[0], &one);
		else if (typeid(T) == typeid(cxfl))   cgemv_ (&trans, &m, &n, &alpha, &A[0], &ah, &x[0], &one, &beta, &y[0], &one);
		else if (typeid(T) == typeid(cxdb))   zgemv_ (&trans, &m, &n, &alpha, &A[0], &ah, &x[0], &one, &beta, &y[0], &one);

		return y;
		
	}
	
};


#endif // __LAPACK_HPP__
