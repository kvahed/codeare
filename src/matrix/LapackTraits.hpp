#include "Complex.hpp"

#include <stdlib.h>
#include "cblas.h"

static int izero = 0;
static int ione  = 1;

static inline int v3p_netlib_dqrsl_ (double *x, int *ldx, int *n, int *k, double *qraux, double *y, 
									 double *qy, double *qty,  double *b, double *rsd, double *xb, 
									 int *job, int *info) {
	return 0; }

static inline int v3p_netlib_dqrdc_ (double* x, int* ldx, int* n, int* p, int* jptv, double* work, int* job) {
	return 0;
}

extern "C" {

	// Cholesky factorization of a complex Hermitian positive definite matrix
	void cpotrf_ (const char* uplo, const int* n, void* a, const int* lda, 
				  int *info);
	void dpotrf_ (const char* uplo, const int* n, void* a, const int* lda, 
				  int *info);
	void spotrf_ (const char* uplo, const int* n, void* a, const int* lda, 
				  int *info);
	void zpotrf_ (const char* uplo, const int* n, void* a, const int* lda, 
				  int *info);
	
	// Computes an LU factorization of a general M-by-N matrix A
	void cgetrf_ (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	void dgetrf_ (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	void zgetrf_ (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	void sgetrf_ (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	
	// Inverse of a complex Hermitian pos def mat with cpotrf/cpptrf
	void cpotri_ (const char* uplo, int*n, void *a, int* lda, int*info);
	void dpotri_ (const char* uplo, int*n, void *a, int* lda, int*info);
	void zpotri_ (const char* uplo, int*n, void *a, int* lda, int*info);
	void spotri_ (const char* uplo, int*n, void *a, int* lda, int*info);
	
	// Matrix inversion through cholesky decomposition
	void cgetri_ (int *n, void *a, int* lda, int *ipiv, void *work, 
				  int *lwork, int*info);
	void dgetri_ (int *n, void *a, int* lda, int *ipiv, void *work, 
				  int *lwork, int*info);
	void zgetri_ (int *n, void *a, int* lda, int *ipiv, void *work, 
				  int *lwork, int*info);
	void sgetri_ (int *n, void *a, int* lda, int *ipiv, void *work, 
				  int *lwork, int*info);
	
	// Eigen value computations
	void cgeev_  (const char* jobvl, const char* jobvr, const int* n, cxfl *a,
                  const int *lda, cxfl *w, cxfl *vl,
                  const int* ldvl, cxfl *vr, const int* ldvr,
				  cxfl *work, const int *lwork, cxfl *rwork, int *info);
	void zgeev_  (const char* jobvl, const char* jobvr, const int* n, cxdb *a,
	              const int *lda, cxdb *w, cxdb *vl,
	              const int* ldvl, cxdb *vr, const int* ldvr,
				  cxdb *work, const int *lwork, cxdb *rwork, int *info);
	void sgeev_  (const char* jobvl, const char* jobvr, const int* n, float *a,
			      const int *lda, float *wr, float *wi, float *vl, const int* ldvl,
			      float *vr, const int* ldvr, float *work, const int *lwork,
			      int *info);
	void dgeev_  (const char* jobvl, const char* jobvr, const int* n, double *a,
			      const int *lda, double *wr, double *wi, double *vl, const int* ldvl,
			      double *vr, const int* ldvr, double *work, const int *lwork,
			      int *info);

	// Singular value decomposition 
	void cgesdd_ (const char *jobz, int *m, int *n, void *a, int *lda, void *s, 
				  void *u, int *ldu, void *vt, int *ldvt, void *work, int *lwork, 
				  void *rwork, int *iwork, int*info);
	void zgesdd_ (const char *jobz, int *m, int *n, void *a, int *lda, void *s, 
				  void *u, int *ldu, void *vt, int *ldvt, void *work, int *lwork, 
				  void *rwork, int *iwork, int*info);
	void dgesdd_ (const char *jobz, int *m, int *n, void *a, int *lda, void *s, 
				  void *u, int *ldu, void *vt, int *ldvt, void *work, int *lwork, 
				  int *iwork, int*info);
	void sgesdd_ (const char *jobz, int *m, int *n, void *a, int *lda, void *s, 
				  void *u, int *ldu, void *vt, int *ldvt, void *work, int *lwork,
				  int *iwork, int*info);
	
	// Pseudo-inversion 
	void zgelsd_ (int* m, int* n, int* nrhs, const void* a, int* lda, void* b, 
				  int* ldb, void* s, void* rcond, int* rank, void* work, 
				  int* lwork, void* rwork, int* iwork, int* info);
	void cgelsd_ (int* m, int* n, int* nrhs, const void* a, int* lda, void* b, 
				  int* ldb, void* s, void* rcond, int* rank, void* work, 
				  int* lwork, void* rwork, int* iwork, int* info);
	void dgelsd_ (int* m, int* n, int* nrhs, const void* a, int* lda, void* b, 
				  int* ldb, void *s, void* rcond, int* rank, void* work, 
				  int* lwork, void* iwork, int* info);
	void sgelsd_ (int* m, int* n, int* nrhs, const void* a, int* lda, void* b, 
				  int* ldb, void *s, void* rcond, int* rank, void* work, 
				  int* lwork, void* iwork, int* info);

	// Matrix vector multiplication
	void sgemv_  (const char* trans, int* m, int* n, void* alpha, const void *a, 
				  int* lda, const void *x, int* incx, void* beta, void *y, 
				  int* incy);
	void dgemv_  (const char* trans, int* m, int* n, void* alpha, const void *a, 
				  int* lda, const void *x, int* incx, void* beta, void *y, 
				  int* incy);
	void cgemv_  (const char* trans, int* m, int* n, void* alpha, const void *a, 
				  int* lda, const void *x, int* incx, void* beta, void *y, 
				  int* incy);
	void zgemv_  (const char* trans, int* m, int* n, void* alpha, const void *a, 
				  int* lda, const void *x, int* incx, void* beta, void *y, 
				  int* incy);

	// Matrix matrix multiplication
	void sgemm_  (const char *transa, const char *transb, int *m, int *n, int *k,
				  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
				  void *beta, void *c, int *ldc);
	void dgemm_  (const char *transa, const char *transb, int *m, int *n, int *k, 
				  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
				  void *beta, void *c, int *ldc);
	void cgemm_  (const char *transa, const char *transb, int *m, int *n, int *k, 
				  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
				  void *beta, void *c, int *ldc);
	void zgemm_  (const char *transa, const char *transb, int *m, int *n, int *k, 
				  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
				  void *beta, void *c, int *ldc);
	
}


template <class T>
struct LapackTraits {};


template <>
struct LapackTraits<float> {

	typedef float Type;

	inline static void 
	potrf (const char* uplo, const int* n, void* a, const int* lda, int *info) {
		spotrf_ (uplo, n, a, lda, info);
	}

	inline static void 
	getrf (int* m, int*n, void *a, int* lda, int*ipiv, int*info) {
		sgetrf_ (m, n, a, lda, ipiv, info);
	}

	inline static void 
	potri (const char* uplo, int*n, void *a, int* lda, int*info) {
		spotri_ (uplo, n, a, lda, info);
	}
	
	inline static void
	getri (int *n, void *a, int* lda, int *ipiv, void *work, int *lwork, 
		   int*info) {
		sgetri_ (n, a, lda, ipiv, work, lwork, info);
	}

	inline static void 
	gemv (const char* trans, int* m, int* n, void* alpha, const void *a, 
		  int* lda, const void *x, int* incx, void* beta, void *y, int* incy) {
		sgemv_ (trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
	}

	inline static void 
	gemm (const char *transa, const char *transb, int  *m, int   *n, int *k, 
		  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
		  void *beta, void *c, int *ldc) {
		sgemm_ (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
	}

    inline static Type 
	nrm2 (const int N, const Type *X, const int incX) {
		return cblas_snrm2 (N, X, incX);
	}
	
	inline static Type
	dot  (const int N, const Type *X, const int incX, const Type *Y, 
		  const int incY, void* res) {
		return cblas_sdot (N, X, incX, Y, incY);	
	}

	inline static Type 
	dotc (const int N, const Type *X, const int incX, const Type *Y, 
		  const int incY, void* res) {
		return cblas_sdot (N, X, incX, Y, incY);	
	}

	inline static void
	geev (const char *jvl, const char *jvr, const int n, Type *a, const int *lda,
		  std::complex<Type> *w, Type *vl, const int *ldvl, Type *vr, const int *ldvr,
		  Type *work, const int *lwork, void* rwork, int *info) {

		int err = 0;
		Type* dwr = (Type*) malloc (n * sizeof(Type));
		Type* dwi = (Type*) malloc (n * sizeof(Type));
		printf ("%d\n",err++);

		Type* dw  = (Type*) w;
		//sgeev_ (jvl, jvr, &n, a, lda, dwr, dwi, vl, ldvl, vr, ldvr, work,
		//		lwork, info);
		printf ("%d\n",err++);

		for (size_t i = 0; i < n; i++) {
			dw[2*i  ] = dwr[i];
			dw[2*i+1] = dwi[i];
		}
		printf ("%d\n",err++);
		free (dwr);
		printf ("%d\n",err++);

		free (dwi);
		printf ("%d\n",err++);

	}

	inline static void 
	gelsd (int* m, int* n, int* nrhs, const void* a, int* lda, void* b, 
		   int* ldb, void* s, Type rcond, int* rank, void* work, 
		   int* lwork, void* rwork, int* iwork, int* info) {
		sgelsd_ (m, n, nrhs, a, lda, b, ldb, s, &rcond, rank, work, lwork, iwork,
				 info);
	}

	inline static void 
	gesdd (const char *jobz, int*m, int *n, void *a, int *lda, void *s, void*u, 
		   int*ldu, void *vt, int *ldvt, void *work, int*lwork, void *rwork, 
		   int *iwork, int*info) {
		sgesdd_ (jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, 
				 info);
	}
		
};


template <>
struct LapackTraits<double> {

	typedef double Type;

	inline static void 
	potrf (const char* uplo, const int* n, void* a, const int* lda, int *info) {
		dpotrf_ (uplo, n, a, lda, info);
	}

	inline static void 
	getrf (int* m, int*n, void *a, int* lda, int*ipiv, int*info) {
		dgetrf_ (m, n, a, lda, ipiv, info);
	}

	inline static void 
	potri (const char* uplo, int*n, void *a, int* lda, int*info) {
		dpotri_ (uplo, n, a, lda, info);
	}
	
	inline static void
	getri (int *n, void *a, int* lda, int *ipiv, void *work, int *lwork, 
		   int*info) {
		dgetri_ (n, a, lda, ipiv, work, lwork, info);
	}

	inline static void 
	gemv (const char* trans, int* m, int* n, void* alpha, const void *a, 
		  int* lda, const void *x, int* incx, void* beta, void *y, int* incy) {
		dgemv_ (trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
	}

	inline static void 
	gemm (const char *transa, const char *transb, int  *m, int   *n, int *k, 
		  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
		  void *beta, void *c, int *ldc) {
		dgemm_ (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
	}

    inline static Type 
	nrm2 (const int N, const Type *X, const int incX) {
		return cblas_dnrm2 (N, X, incX);
	}
	
	inline static Type 
	dot  (const int N, const Type *X, const int incX, const Type *Y, 
		  const int incY, void* res) {
		return cblas_ddot (N, X, incX, Y, incY);	
	}

	inline static Type
	dotc (const int N, const Type *X, const int incX, const Type *Y, 
		  const int incY, void* res) {
		return cblas_ddot (N, X, incX, Y, incY);	
	}

	inline static void
	geev (const char *jvl, const char *jvr, const int n, Type *a, const int *lda,
		  std::complex<Type> *w, Type *vl, const int *ldvl, Type *vr, const int *ldvr,
		  Type *work, const int *lwork, void* rwork, int *info) {

		Type* dwr = (Type*) malloc (n * sizeof(Type));
		Type* dwi = (Type*) malloc (n * sizeof(Type));
		Type* dw  = (Type*) w;

		dgeev_ (jvl, jvr, &n, a, lda, dwr, dwi, vl, ldvl, vr, ldvr, work, 
				lwork, info);
				
		for (size_t i = 0; i < n; i++) {
			dw[2*i  ] = dwr[i];
			dw[2*i+1] = dwi[i];
		}

		free (dwr);
		free (dwi);

	}

	inline static void 
	gelsd (int* m, int* n, int* nrhs, const void* a, int* lda, void* b, 
		   int* ldb, void* s, Type rcond, int* rank, void* work, 
		   int* lwork, void* rwork, int* iwork, int* info) {
		dgelsd_ (m, n, nrhs, a, lda, b, ldb, s, &rcond, rank, work, lwork, iwork, 
				 info);
	}

	inline static void 
	gesdd (const char *jobz, int*m, int *n, void *a, int *lda, void *s, void*u, 
		   int*ldu, void *vt, int *ldvt, void *work, int*lwork, void *rwork, 
		   int *iwork, int*info) {
		dgesdd_ (jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, 
				 info);
	}
		
};



template <>
struct LapackTraits<cxfl> {

	typedef cxfl Type;

	inline static void 
	potrf (const char* uplo, const int* n, void* a, const int* lda, int *info) {
		cpotrf_ (uplo, n, a, lda, info);
	}

	inline static void 
	getrf (int* m, int*n, void *a, int* lda, int*ipiv, int*info) {
		cgetrf_ (m, n, a, lda, ipiv, info);
	}

	inline static void 
	potri (const char* uplo, int*n, void *a, int* lda, int*info) {
		cpotri_ (uplo, n, a, lda, info);
	}
	
	inline static void
	getri (int *n, void *a, int* lda, int *ipiv, void *work, int *lwork, 
		   int*info) {
		cgetri_ (n, a, lda, ipiv, work, lwork, info);
	}

	inline static void 
	gemv (const char* trans, int* m, int* n, void* alpha, const void *a, 
		  int* lda, const void *x, int* incx, void* beta, void *y, int* incy) {
		cgemv_ (trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
	}

	inline static void 
	gemm (const char *transa, const char *transb, int  *m, int   *n, int *k, 
		  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
		  void *beta, void *c, int *ldc) {
		cgemm_ (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
	}

    inline static Type
	nrm2 (const int N, const void *X, const int incX) {
		return cblas_scnrm2 (N, X, incX);
	}
	
	inline static void 
	dot  (const int N, const void *X, const int incX, const void *Y, 
		  const int incY, void* res) {
		cblas_cdotu_sub (N, X, incX, Y, incY, res);	
	}

	inline static void 
	dotc (const int N, const void *X, const int incX, const void *Y, 
		  const int incY, void* res) {
		cblas_cdotc_sub (N, X, incX, Y, incY, res);	
	}

	inline static void 
	geev (const char *jvl, const char *jvr, const int& n, Type *a, const int *lda,
		  Type *w, Type *vl, const int *ldvl, Type *vr, const int *ldvr, Type *work,
		  const int *lwork, Type* rwork, int *info) {
		cgeev_  (jvl, jvr, &n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, 
				 rwork, info);
	}

	inline static void 
	gelsd (int* m, int* n, int* nrhs, const void* a, int* lda, void* b, 
		   int* ldb, void* s, double rcond, int* rank, void* work, int* lwork, 
		   void* rwork, int* iwork, int* info) {
		float frcond = float(rcond);
		cgelsd_ (m, n, nrhs, a, lda, b, ldb, s, &frcond, rank, work, lwork, 
				 rwork, iwork, info);
	}

	inline static void 
	gesdd (const char *jobz, int*m, int *n, void *a, int *lda, void *s, void *u, 
		   int*ldu, void *vt, int *ldvt, void *work, int*lwork, void *rwork, 
		   int *iwork, int*info) {
		cgesdd_ (jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, 
				 iwork, info);
	}

};



template <>
struct LapackTraits<cxdb> {

	typedef cxdb Type;

	inline static void 
	potrf (const char* uplo, const int* n, void* a, const int* lda, int *info) {
		zpotrf_ (uplo, n, a, lda, info);
	}

	inline static void 
	getrf (int* m, int*n, void *a, int* lda, int*ipiv, int*info) {
		zgetrf_ (m, n, a, lda, ipiv, info);
	}

	inline static void 
	potri (const char* uplo, int*n, void *a, int* lda, int*info) {
		zpotri_ (uplo, n, a, lda, info);
	}
	
	inline static void
	getri (int *n, void *a, int* lda, int *ipiv, void *work, int *lwork, 
		   int*info) {
		zgetri_ (n, a, lda, ipiv, work, lwork, info);
	}

	inline static void 
	gemv (const char* trans, int* m, int* n, void* alpha, const void *a, 
		  int* lda, const void *x, int* incx, void* beta, void *y, int* incy) {
		zgemv_ (trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
	}

	inline static void 
	gemm (const char *transa, const char *transb, int  *m, int   *n, int *k, 
		  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
		  void *beta, void *c, int *ldc) {
		zgemm_ (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
	}

    inline static Type 
	nrm2 (const int N, const void *X, const int incX) {
		return cblas_dznrm2 (N, X, incX);
	}
	
	inline static void 
	dot  (const int N, const void *X, const int incX, const void *Y, 
		  const int incY, void* res) {
		cblas_zdotu_sub (N, X, incX, Y, incY, res);	
	}

	inline static void 
	dotc (const int N, const void *X, const int incX, const void *Y, 
		  const int incY, void* res) {
		cblas_zdotc_sub (N, X, incX, Y, incY, res);	
	}

	inline static void
	geev (const char *jvl, const char *jvr, const int n, Type *a, const int *lda,
		  Type *w, Type *vl, const int *ldvl, Type *vr, const int *ldvr, Type *work,
		  int *lwork, Type* rwork, int *info) {
		zgeev_  (jvl, jvr, &n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, 
				 rwork, info);
	}

	inline static void 
	gelsd (int* m, int* n, int* nrhs, const void* a, int* lda, void* b, 
		   int* ldb, void* s, double rcond, int* rank, void* work, int* lwork, 
		   void* rwork, int* iwork, int* info) {
		zgelsd_ (m, n, nrhs, a, lda, b, ldb, s, &rcond, rank, work, lwork, rwork, 
				 iwork, info);
	}

	inline static void 
	gesdd (const char *jobz, int*m, int *n, void *a, int *lda, void *s, void *u, 
		   int*ldu, void *vt, int *ldvt, void *work, int*lwork, void *rwork, 
		   int *iwork, int*info) {
		zgesdd_ (jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, 
				 iwork, info);
	}


};


