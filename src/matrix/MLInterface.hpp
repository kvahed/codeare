#ifdef USE_ACML

  // Declare wrappers for ACML versions of BLAS/LAPACK routines
  #undef  UNDERSCORE
  #define NOUNDERSCORE
  #define F77name
  #include "acml.h"

#else

  // NBN: Win32 ATLAS, MKL
  #if defined(WIN32) && !defined(__CYGWIN__)
    #ifndef UNDERSCORE
    #define UNDERSCORE
    #endif
  #elif defined(USE_MKL)
    #ifndef NOUNDERSCORE
    #define NOUNDERSCORE
    #endif
  #endif

  #if defined(MRVISTA)
    #ifndef NOUNDERSCORE
    #define NOUNDERSCORE
    #endif
  #endif

  #if defined(__APPLE__)
    #ifndef UNDERSCORE
    #define UNDERSCORE
    #endif
  #endif


  #if defined(UNDERSCORE)
    #define F77name(X,UX) X ## _
  #elif defined(PREUNDERSCORE)
    #define F77name(X,UX) _ ## X ## _
  #elif defined(NOUNDERSCORE)
    #define F77name(X,UX) X
  #endif

extern "C" {

	// Cholesky factorization of a complex Hermitian positive definite matrix
	void F77name(spotrf,CPOTRF) (const char* uplo, const int* n, void* a, const int* lda, 
				  int *info);
	void F77name(dpotrf,DPOTRF) (const char* uplo, const int* n, void* a, const int* lda, 
				  int *info);
	void F77name(cpotrf,CPOTRF) (const char* uplo, const int* n, void* a, const int* lda, 
				  int *info);
	void F77name(zpotrf,ZPOTRF) (const char* uplo, const int* n, void* a, const int* lda, 
				  int *info);
	
	// Computes an LU factorization of a general M-by-N matrix A
	void F77name(sgetrf,SGETRF) (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	void F77name(dgetrf,DGETRF) (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	void F77name(cgetrf,CGETRF) (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	void F77name(zgetrf,ZGETRF) (int* m, int*n, void *a, int* lda, int*ipiv, int*info);
	
	// Inverse of a complex Hermitian pos def mat with cpotrf/cpptrf
	void F77name(spotri,SPOTRI) (const char* uplo, int*n, void *a, int* lda, int*info);
	void F77name(dpotri,DPOTRI) (const char* uplo, int*n, void *a, int* lda, int*info);
	void F77name(cpotri,CPOTRI) (const char* uplo, int*n, void *a, int* lda, int*info);
	void F77name(zpotri,ZPOTRI) (const char* uplo, int*n, void *a, int* lda, int*info);
	
	// Matrix inversion through cholesky decomposition
	void F77name(sgetri,SGETRI) (int *n, void *a, int* lda, int *ipiv, void *work, 
				  int *lwork, int*info);
	void F77name(dgetri,DGETRI) (int *n, void *a, int* lda, int *ipiv, void *work, 
				  int *lwork, int*info);
	void F77name(cgetri,CGETRI) (int *n, void *a, int* lda, int *ipiv, void *work, 
				  int *lwork, int*info);
	void F77name(zgetri,ZGETRI) (int *n, void *a, int* lda, int *ipiv, void *work, 
				  int *lwork, int*info);
	
	// Eigen value computations
	void F77name(sgeev,SGEEV)  (const char *jvl, const char *jvr, int *n, const void *a, 
				  int *lda, void *wr, void *wi, void *vl, int *ldvl, void *vr, 
				  int *ldvr, void *work, int *lwork, int *info);
	void F77name(dgeev,DGEEV)  (const char *jvl, const char *jvr, int *n, const void *a, 
				  int *lda, void *wr, void *wi, void *vl, int *ldvl, void *vr, 
				  int *ldvr, void *work, int *lwork, int *info);
	void F77name(cgeev,CGEEV)  (const char *jvl, const char *jvr, int *n, const void *a, 
				  int *lda, void *w,  void *vl, int *ldvl, void *vr, int *ldvr, 
				  void *work, int *lwork, void *rwork, int *info);
	void F77name(zgeev,ZGEEV)  (const char *jvl, const char *jvr, int *n, const void *a, 
				  int *lda, void *w,  void *vl, int *ldvl, void *vr, int *ldvr, 
				  void *work, int *lwork, void *rwork, int *info);
	
	// Singular value decomposition 
	void F77name(sgesdd,SGESDD) (const char *jobz, int*m, int *n, void *a, int *lda, void *s, 
				  void*u, int*ldu, void *vt, int *ldvt, void *work, int*lwork, 
				  int *iwork, int*info);
	void F77name(dgesdd,DGESDD) (const char *jobz, int*m, int *n, void *a, int *lda, void *s, 
				  void*u, int*ldu, void *vt, int *ldvt, void *work, int*lwork, 
				  int *iwork, int*info);
	void F77name(cgesdd,CGESDD) (const char *jobz, int*m, int *n, void *a, int *lda, void *s, 
				  void*u, int*ldu, void *vt, int *ldvt, void *work, int*lwork, 
				  void *rwork, int *iwork, int*info);
	void F77name(zgesdd,ZGESDD) (const char *jobz, int*m, int *n, void *a, int *lda, void *s, 
				  void*u, int*ldu, void *vt, int *ldvt, void *work, int*lwork, 
				  void *rwork, int *iwork, int*info);
	
	// Pseudo-inversion 
	void F77name(sgels,SGELS) (const char* trans, const int* m, const int* n, const int* nrhs,
			float* a, const int* lda, float* b, const int* ldb, float* work, const int* lwork,
			int* info);
	void F77name(dgels,DGELS) (const char* trans, const int* m, const int* n, const int* nrhs,
			double* a, const int* lda, double* b, const int* ldb, double* work, const int* lwork,
			int* info);
	void F77name(cgels,CGELS) (const char* trans, const int* m, const int* n, const int* nrhs,
			cxfl* a, const int* lda, cxfl* b, const int* ldb, cxfl* work, const int* lwork,
			int* info);
	void F77name(zgels,ZGELS) (const char* trans, const int* m, const int* n, const int* nrhs,
			cxdb* a, const int* lda, cxdb* b, const int* ldb, cxdb* work, const int* lwork,
			int* info);

	// Pseudo-inversion
		void F77name(sgelsd,SGELSD) (int* m, int* n, int* nrhs, float* a, int* lda, void* b,
					  int* ldb, void *s, void* rcond, int* rank, void* work,
					  int* lwork, void* iwork, int* info);
		void F77name(dgelsd,DGELSD) (int* m, int* n, int* nrhs, double* a, int* lda, void* b,
					  int* ldb, void *s, void* rcond, int* rank, void* work,
					  int* lwork, void* iwork, int* info);
		void F77name(cgelsd,CGELSD) (int* m, int* n, int* nrhs, cxfl* a, int* lda, void* b,
					  int* ldb, void* s, void* rcond, int* rank, void* work,
					  int* lwork, void* rwork, int* iwork, int* info);
		void F77name(zgelsd,ZGELSD) (int* m, int* n, int* nrhs, cxdb* a, int* lda, void* b,
					  int* ldb, void* s, void* rcond, int* rank, void* work,
					  int* lwork, void* rwork, int* iwork, int* info);

		// Matrix vector multiplication
	void F77name(sgemv,SGEMV)  (const char* trans, int* m, int* n, void* alpha, const void *a, 
				  int* lda, const void *x, int* incx, void* beta, void *y, 
				  int* incy);
	void F77name(dgemv,DGEMV)  (const char* trans, int* m, int* n, void* alpha, const void *a, 
				  int* lda, const void *x, int* incx, void* beta, void *y, 
				  int* incy);
	void F77name(cgemv,CGEMV)  (const char* trans, int* m, int* n, void* alpha, const void *a, 
				  int* lda, const void *x, int* incx, void* beta, void *y, 
				  int* incy);
	void F77name(zgemv,ZGEMV)  (const char* trans, int* m, int* n, void* alpha, const void *a, 
				  int* lda, const void *x, int* incx, void* beta, void *y, 
				  int* incy);

	// Matrix matrix multiplication
	void F77name(sgemm,SGEMM)  (const char *transa, const char *transb, int *m, int *n, int *k,
				  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
				  void *beta, void *c, int *ldc);
	void F77name(dgemm,DGEMM)  (const char *transa, const char *transb, int *m, int *n, int *k, 
				  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
				  void *beta, void *c, int *ldc);
	void F77name(cgemm,CGEMM)  (const char *transa, const char *transb, int *m, int *n, int *k, 
				  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
				  void *beta, void *c, int *ldc);
	void F77name(zgemm,ZGEMM) (const char *transa, const char *transb, int *m, int *n, int *k, 
				  void *alpha, const void *a, int *lda, const void *b, int *ldb, 
				  void *beta, void *c, int *ldc);

}	

#endif

#ifdef _MSC_VER
extern "C" long _ftol( double ); //defined by VC6 C libs
extern "C" long _ftol2( double dblSource ) { return _ftol( dblSource ); }
#endif

#define SGEEV  F77name(sgeev,SGEEV)
#define SGETRF F77name(sgetrf,SGETRF)
#define SGETRI F77name(sgetri,SGETRI) 
#define SPOTRF F77name(spotrf,SPOTRF)
#define SPOTRI F77name(spotri,SPOTRI)
#define SGELS  F77name(sgels,SGELS)
#define SGELSD F77name(sgelsd,SGELSD) 
#define SGESDD F77name(sgesdd,SGESDD)
#define SGEMM  F77name(sgemm,SGEMM) 
#define SGEMV  F77name(sgemv,SGEMV) 

#define DGEEV  F77name(dgeev,DGEEV)
#define DGETRF F77name(dgetrf,DGETRF)
#define DGETRI F77name(dgetri,DGETRI) 
#define DPOTRF F77name(dpotrf,DPOTRF)
#define DPOTRI F77name(dpotri,DPOTRI)
#define DGELS  F77name(dgels,DGELS)
#define DGELSD F77name(dgelsd,DGELSD) 
#define DGESDD F77name(dgesdd,DGESDD)
#define DGEMM  F77name(dgemm,DGEMM) 
#define DGEMV  F77name(dgemv,DGEMV) 

#define CGEEV  F77name(cgeev,CGEEV)
#define CGETRF F77name(cgetrf,CGETRF)
#define CGETRI F77name(cgetri,CGETRI) 
#define CPOTRF F77name(cpotrf,CPOTRF)
#define CPOTRI F77name(cpotri,CPOTRI)
#define CGELS  F77name(cgels,CGELS)
#define CGELSD F77name(cgelsd,CGELSD) 
#define CGESDD F77name(cgesdd,CGESDD)
#define CGEMM  F77name(cgemm,CGEMM) 
#define CGEMV  F77name(cgemv,CGEMV) 

#define ZGEEV  F77name(zgeev,ZGEEV)
#define ZGETRF F77name(zgetrf,ZGETRF)
#define ZGETRI F77name(zgetri,ZGETRI) 
#define ZPOTRF F77name(zpotrf,ZPOTRF)
#define ZPOTRI F77name(zpotri,ZPOTRI)
#define ZGELS  F77name(zgels,ZGELS)
#define ZGELSD F77name(zgelsd,ZGELSD) 
#define ZGESDD F77name(zgesdd,ZGESDD)
#define ZGEMM  F77name(zgemm,ZGEMM) 
#define ZGEMV  F77name(zgemv,ZGEMV) 
