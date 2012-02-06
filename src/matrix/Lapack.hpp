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

}

//#include "Matrix.hpp"

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
		
		int    N     = m.Height();

		int    lda   =  N;
		int    ldvl  =  (jobvl) ? N : 1;
		int    ldvr  =  (jobvr) ? N : 1;
		int    info  =  0;
		int    lwork = -1;
		
		T*     w;
		T*     wi;

		if (typeid(T) == typeid(float) || typeid(T) == typeid(double)) {
			w  = (T*)     malloc (  N * sizeof(T));
			wi = (T*)     malloc (  N * sizeof(T));
		}

		T*     work  = (T*)     malloc (  1 * sizeof(T));
		float* rwork = (float*) malloc (2*N * sizeof(float));
		
		// Workspace query
		if (typeid(T) == typeid(cxfl))
			cgeev_ (&jobvl, &jobvr, &N, &m[0], &lda, &ev[0], &lv[0], &ldvl, &rv[0], &ldvr, work, &lwork, rwork, &info);
		else if (typeid(T) == typeid(cxdb))
			zgeev_ (&jobvl, &jobvr, &N, &m[0], &lda, &ev[0], &lv[0], &ldvl, &rv[0], &ldvr, work, &lwork, rwork, &info);
		else if (typeid(T) == typeid(double))
			dgeev_ (&jobvl, &jobvr, &N, &m[0], &lda,  w, wi, &lv[0], &ldvl, &rv[0], &ldvr, work, &lwork,        &info);
		else if (typeid(T) == typeid(float))
			sgeev_ (&jobvl, &jobvr, &N, &m[0], &lda,  w, wi, &lv[0], &ldvl, &rv[0], &ldvr, work, &lwork,        &info);
		
		// Intialise work space
		lwork = (int) (cxfl(work[0]).real());
		free (work); work = (T*) malloc (lwork*sizeof(T));
		
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
		}
		free (work);
		free (rwork);
		
		return info;
		
	}
	

	template<class T, class S> static int 
	SVD (Matrix<T>& m, Matrix<S>& s, Matrix<T>& u, Matrix<T>& v, const char jobz = 'N') {
		
		// SVD only defined on 2D data
		if (!m.Is2D())
			return -2;
		
		bool   cm    = true;
		
		int    M     = m.Height();
		int    N     = m.Width();
		int    lwork = -1;
		int    info  = 0;
		int    lda   = M;
		int    ldu   = MAX(M,N);
		int    ldvt  = MIN(M,N);
		
		T*     work  =     (T*) malloc (sizeof(T));
		float* rwork = (float*) malloc (sizeof(float));
		int*   iwork =   (int*) malloc (8 * MIN(M,N) * sizeof(int));
		
		// Only needed for complex data
		if (typeid(T) == typeid(cxfl)) {
			free (rwork);
			rwork = (float*) malloc (MIN(M,N) * MAX(5*MIN(M,N)+7,2*MAX(M,N)+2*MIN(M,N)+1) * sizeof(float));
		}
		
		// Workspace query
		if (typeid(T) == typeid(cxfl))
			cgesdd_ (&jobz, &M, &N, &m[0], &lda, (float*)&s[0], &u[0], &ldu, &v[0], &ldvt, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(double))
			dgesdd_ (&jobz, &M, &N, &m[0], &lda,         &s[0], &u[0], &ldu, &v[0], &ldvt, work, &lwork,        iwork, &info);
		
		// Resize work according to ws query
		lwork = (int) cxfl (work[0]).real();
		free(work); work = (T*) malloc (lwork * sizeof(T));
		
		//SVD
		if (typeid(T) == typeid(cxfl))
			cgesdd_ (&jobz, &M, &N, &m[0], &lda, (float*)&s[0], &u[0], &ldu, &v[0], &ldvt, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(double))
			dgesdd_ (&jobz, &M, &N, &m[0], &lda,         &s[0], &u[0], &ldu, &v[0], &ldvt, work, &lwork,        iwork, &info);

		v = v.tr();

		// Clean up
		free (work);
		free (rwork);
		free (iwork);
		
		return info;
		
	} 
	

	template <class T> static Matrix<T> 
	Inv (Matrix<T>& m) {
		
		// 2D 
		if (!m.Is2D()) {
			printf ("Inv Error: Parameter m must be 2D");
			return -2;
		}

		// Square matrix
		if (m.Width() != m.Height()){
			printf ("Inv Error: Parameter m must be square");
			return -3;
		}
		
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
		
		int  M      =  m.Height();
		int  N      =  m.Width();
		int  nrhs   =  M;
		int  lda    =  M;
		int  ldb    =  MAX(M,N);
		int  lwork  = -1; 
		int  rank   =  0;
		int  info   =  0;
		
		T*   work   = (T*)   malloc (         sizeof(T)); 
		int* iwork  = (int*) malloc (         sizeof(int)); 
		T*   rwork  = (T*)   malloc (         sizeof(float));    
		T*   s      = (T*)   malloc (MIN(M,N)*sizeof(float));    
		
		Matrix<T> b =  Matrix<T>::Id(ldb);
		
		float     rcond  = -1.0;
		
		if (typeid(T) == typeid(cxfl))
			cgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &rcond, &rank, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(cxdb))
			zgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &rcond, &rank, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(double))
			dgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &rcond, &rank, work, &lwork,        iwork, &info);
		else if (typeid(T) == typeid(float))
			sgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &rcond, &rank, work, &lwork,        iwork, &info);
		
		if (typeid(T) == typeid(cxfl))
			rwork = (float*) realloc (rwork, (int)rwork[0]*sizeof(float));
		
		iwork = (int*) realloc (iwork, iwork[0]*sizeof(int));

		if (typeid (T) == typeid(cxfl) || typeid (T) == typeid(float))
			lwork = (int) ((float*)work)[0];
		else if (typeid (T) == typeid(cxdb) || typeid (T) == typeid(double))
			lwork = (int) ((double*)work)[0];

		work  = (T*) realloc (work, lwork * sizeof(T));
		
		if (typeid(T) == typeid(cxfl))
			cgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &rcond, &rank, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(cxdb))
			zgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &rcond, &rank, work, &lwork, rwork, iwork, &info);
		else if (typeid(T) == typeid(double))
			dgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &rcond, &rank, work, &lwork,        iwork, &info);
		else if (typeid(T) == typeid(float))
			sgelsd_ (&M, &N, &nrhs, &m[0], &lda, &b[0], &ldb, s, &rcond, &rank, work, &lwork,        iwork, &info);
		
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
		
		//dpotrf_ ( &uplo, &tmp.Width(), &tmp[0], &tmp.Height(), &info);
		
		return res;
		
	}
	
};
