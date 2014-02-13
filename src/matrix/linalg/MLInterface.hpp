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
  #elif defined(MRVISTA)
    #ifndef NOUNDERSCORE
    #define NOUNDERSCORE
    #endif
  #else
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
	void F77name(spotrf,CPOTRF) (const char* uplo, const int* n,  float* a, const int* lda, int *info);
	void F77name(dpotrf,DPOTRF) (const char* uplo, const int* n, double* a, const int* lda, int *info);
	void F77name(cpotrf,CPOTRF) (const char* uplo, const int* n,   cxfl* a, const int* lda, int *info);
	void F77name(zpotrf,ZPOTRF) (const char* uplo, const int* n,   cxdb* a, const int* lda, int *info);
	
	// Computes an LU factorization of a general M-by-N matrix A
	void F77name(sgetrf,SGETRF) (const int* m, const int *n,  float *a, const int* lda, int *ipiv, int *info);
	void F77name(dgetrf,DGETRF) (const int* m, const int *n, double *a, const int* lda, int *ipiv, int *info);
	void F77name(cgetrf,CGETRF) (const int* m, const int *n,   cxfl *a, const int* lda, int *ipiv, int *info);
	void F77name(zgetrf,ZGETRF) (const int* m, const int *n,   cxdb *a, const int* lda, int *ipiv, int *info);
	
	// Inverse of a complex Hermitian pos def mat with cpotrf/cpptrf
	void F77name(spotri,SPOTRI) (const char* uplo, int*n, void *a, int* lda, int*info);
	void F77name(dpotri,DPOTRI) (const char* uplo, int*n, void *a, int* lda, int*info);
	void F77name(cpotri,CPOTRI) (const char* uplo, int*n, void *a, int* lda, int*info);
	void F77name(zpotri,ZPOTRI) (const char* uplo, int*n, void *a, int* lda, int*info);
	
	// Inverse of a complex Hermitian pos def mat with cpotrf/cpptrf
#ifdef __APPLE__
    double F77name(sdot, SDOT)  (const int* n, const float* x, const int* incx, const float* y, const int* incy);
    void   F77name(cdotc,CDOTC) (cxfl*, const int* n, const cxfl* x, const int* incx, const cxfl* y, const int* incy);
    void   F77name(cdotu,CDOTU) (cxfl*, const int* n, const cxfl* x, const int* incx, const cxfl* y, const int* incy);
    void   F77name(zdotc,ZDOTC) (cxdb*, const int* n, const cxdb* x, const int* incx, const cxdb* y, const int* incy);
    void   F77name(zdotu,ZDOTU) (cxdb*, const int* n, const cxdb* x, const int* incx, const cxdb* y, const int* incy);
#else
    float  F77name(sdot, SDOT)  (const int* n, const float* x, const int* incx, const float* y, const int* incy);
    cxfl   F77name(cdotc,CDOTC) (const int* n, const cxfl* x, const int* incx, const cxfl* y, const int* incy);
    cxfl   F77name(cdotu,CDOTU) (const int* n, const cxfl* x, const int* incx, const cxfl* y, const int* incy);
    cxdb   F77name(zdotc,ZDOTC) (const int* n, const cxdb* x, const int* incx, const cxdb* y, const int* incy);
    cxdb   F77name(zdotu,ZDOTU) (const int* n, const cxdb* x, const int* incx, const cxdb* y, const int* incy);
#endif
    double F77name(ddot, DDOT)  (const int* n, const double* x, const int* incx, const double* y, const int* incy);

#ifdef __APPLE__
    double F77name(snrm2,SNRM2) (const int* n, const float* x, const int* incx);
    double F77name(scnrm2,SCNRM2) (const int* n, const cxfl* x, const int* incx);
#else
    float F77name(snrm2,SNRM2) (const int* n, const float* x, const int* incx);
    float F77name(scnrm2,SCNRM2) (const int* n, const cxfl* x, const int* incx);
#endif
    double F77name(dnrm2,DNRM2) (const int* n, const double* x, const int* incx);
    double F77name(dznrm2,DZNRM2) (const int* n, const cxdb* x, const int* incx);

	// Matrix inversion through cholesky decomposition
	void F77name(sgetri,SGETRI) (const int *n,  float *a, const int* lda, int *ipiv,  float *work, const int *lwork,
                                 int* info);
	void F77name(dgetri,DGETRI) (const int *n, double *a, const int* lda, int *ipiv, double *work, const int *lwork,
                                 int* info);
	void F77name(cgetri,CGETRI) (const int *n,   cxfl *a, const int* lda, int *ipiv,   cxfl *work, const int *lwork,
                                 int* info);
	void F77name(zgetri,ZGETRI) (const int *n,   cxdb *a, const int* lda, int *ipiv,   cxdb *work, const int *lwork,
                                 int* info);
	
	// Eigen value computations
	void F77name(sgeev,SGEEV)  (const char* jvl, const char* jvr, const int* n,  float *a, const int* lda,  float *wr,
                                 float *wi,  float *vl, const int* ldvl,  float *vr, const int* ldvr,  float *work,
                                const int* lwork, int* info);
	void F77name(dgeev,DGEEV)  (const char* jvl, const char* jvr, const int* n, double *a, const int* lda, double *wr,
                                double *wi, double *vl, const int* ldvl, double *vr, const int* ldvr, double *work,
                                const int* lwork, int* info);
	void F77name(cgeev,CGEEV)  (const char* jvl, const char* jvr, const int *n,   cxfl *a, const int *lda,
                                  cxfl *w,    cxfl *vl, const int *ldvl,   cxfl *vr, const int *ldvr,   cxfl *work,
                                const int *lwork,  float *rwork, int *info);
	void F77name(zgeev,ZGEEV)  (const char* jvl, const char* jvr, const int *n,   cxdb *a, const int *lda,
                                  cxdb *w,   cxdb *vl, const int *ldvl,   cxdb *vr, const int *ldvr,    cxdb *work,
                                const int *lwork, double *rwork, int *info);
	
	// Singular value decomposition 
	void F77name(sgesdd,SGESDD) (const char* jobz, const int *m, const int *n,  float *a, const int *lda,  float *s,
                                  float* u, const int *ldu,  float *vt, const int *ldvt,  float *work, const int *lwork,
                                 int *iwork, int *info);
	void F77name(dgesdd,DGESDD) (const char* jobz, const int *m, const int *n, double *a, const int *lda, double *s,
                                 double* u, const int *ldu, double *vt, const int *ldvt, double *work, const int *lwork,
                                 int *iwork, int *info);
    void F77name(cgesdd,CGESDD) (const char *jobz, const int *m, const int *n,   cxfl *a, const int *lda,  float *s,
                                   cxfl* u, const int* ldu,   cxfl *vt, const int *ldvt,   cxfl *work, const int* lwork,
                                  float *rwork, int *iwork, int* info);
	void F77name(zgesdd,ZGESDD) (const char *jobz, const int *m, const int *n,   cxdb *a, const int *lda, double *s,
                                   cxdb* u, const int* ldu,   cxdb *vt, const int *ldvt,   cxdb *work, const int* lwork,
                                 double *rwork, int *iwork, int* info);
    
	// Pseudo-inversion 
	void F77name(sgels,SGELS) (const char* trans, const int* m, const int* n, const int* nrhs,  float* a, const int* lda,
                                float* b, const int* ldb,  float* work, const int* lwork, int* info);
	void F77name(dgels,DGELS) (const char* trans, const int* m, const int* n, const int* nrhs, double* a, const int* lda,
                               double* b, const int* ldb, double* work, const int* lwork, int* info);
	void F77name(cgels,CGELS) (const char* trans, const int* m, const int* n, const int* nrhs,   cxfl* a, const int* lda,
                                 cxfl* b, const int* ldb,   cxfl* work, const int* lwork, int* info);
	void F77name(zgels,ZGELS) (const char* trans, const int* m, const int* n, const int* nrhs,   cxdb* a, const int* lda,
                                 cxdb* b, const int* ldb,   cxdb* work, const int* lwork, int* info);
    
	// Matrix vector multiplication
	void F77name(sgemv,SGEMV)  (const char* trans, const int* m, const int* n, const  float* alpha, const  float *a,
                                const int* lda, const  float *x, const int* incx, const float* beta,   float *y,
                                const int* incy);
	void F77name(dgemv,DGEMV)  (const char* trans, const int* m, const int* n, const double* alpha, const double *a,
                                const int* lda, const double *x, const int* incx, const double* beta, double *y,
                                const int* incy);
	void F77name(cgemv,CGEMV)  (const char* trans, const int* m, const int* n, const   cxfl* alpha, const   cxfl *a,
                                const int* lda, const   cxfl *x, const int* incx, const  cxfl* beta,    cxfl *y,
                                const int* incy);
	void F77name(zgemv,ZGEMV)  (const char* trans, const int* m, const int* n, const   cxdb* alpha, const   cxdb *a,
                                const int* lda, const   cxdb *x, const int* incx, const  cxdb* beta,    cxdb *y,
                                const int* incy);

	// Matrix matrix multiplication
	void F77name(sgemm,SGEMM)  (const char *transa, const char *transb, const int *m, const int *n, const int *k,
                                const  float *alpha, const  float *a, const int *lda, const  float *b, const int *ldb,
                                const  float *beta,  float *c, const int *ldc);
	void F77name(dgemm,DGEMM)  (const char *transa, const char *transb, const int *m, const int *n, const int *k,
                                const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
                                const double *beta, double *c, const int *ldc);
	void F77name(cgemm,CGEMM)  (const char *transa, const char *transb, const int *m, const int *n, const int *k,
                                const   cxfl *alpha, const   cxfl *a, const int *lda, const   cxfl *b, const int *ldb,
                                const   cxfl *beta,   cxfl *c, const int *ldc);
	void F77name(zgemm,ZGEMM)  (const char *transa, const char *transb, const int *m, const int *n, const int *k,
                                const   cxdb *alpha, const   cxdb *a, const int *lda, const   cxdb *b, const int *ldb,
                                const   cxdb *beta,   cxdb *c, const int *ldc);
    
}	

#endif

#ifdef _MSC_VER
extern "C" long _ftol( double ); //defined by VC6 C libs
extern "C" long _ftol2( double dblSource ) { return _ftol( dblSource ); }
#endif

#define SDOT   F77name(sdot,SDOT)
#define SNRM2  F77name(snrm2,SNRM2)
#define SGEEV  F77name(sgeev,SGEEV)
#define SGETRF F77name(sgetrf,SGETRF)
#define SGETRI F77name(sgetri,SGETRI) 
#define SPOTRF F77name(spotrf,SPOTRF)
#define SPOTRI F77name(spotri,SPOTRI)
#define SGELS  F77name(sgels,SGELS)
#define SGESDD F77name(sgesdd,SGESDD)
#define SGEMM  F77name(sgemm,SGEMM) 
#define SGEMV  F77name(sgemv,SGEMV)

#define DDOT   F77name(ddot,DDOT)
#define DNRM2  F77name(dnrm2,DNRM2)
#define DGEEV  F77name(dgeev,DGEEV)
#define DGETRF F77name(dgetrf,DGETRF)
#define DGETRI F77name(dgetri,DGETRI) 
#define DPOTRF F77name(dpotrf,DPOTRF)
#define DPOTRI F77name(dpotri,DPOTRI)
#define DGELS  F77name(dgels,DGELS)
#define DGESDD F77name(dgesdd,DGESDD)
#define DGEMM  F77name(dgemm,DGEMM) 
#define DGEMV  F77name(dgemv,DGEMV) 

#define CDOTU  F77name(cdotu,CDOTU)
#define CDOTC  F77name(cdotc,CDOTC)
#define SCNRM2 F77name(scnrm2,SCNRM2)
#define CGEEV  F77name(cgeev,CGEEV)
#define CGETRF F77name(cgetrf,CGETRF)
#define CGETRI F77name(cgetri,CGETRI) 
#define CPOTRF F77name(cpotrf,CPOTRF)
#define CPOTRI F77name(cpotri,CPOTRI)
#define CGELS  F77name(cgels,CGELS)
#define CGESDD F77name(cgesdd,CGESDD)
#define CGEMM  F77name(cgemm,CGEMM) 
#define CGEMV  F77name(cgemv,CGEMV) 

#define ZDOTU  F77name(zdotu,ZDOTU)
#define ZDOTC  F77name(zdotc,ZDOTC)
#define DZNRM2 F77name(dznrm2,DZNRM2)
#define ZGEEV  F77name(zgeev,ZGEEV)
#define ZGETRF F77name(zgetrf,ZGETRF)
#define ZGETRI F77name(zgetri,ZGETRI) 
#define ZPOTRF F77name(zpotrf,ZPOTRF)
#define ZPOTRI F77name(zpotri,ZPOTRI)
#define ZGELS  F77name(zgels,ZGELS)
#define ZGESDD F77name(zgesdd,ZGESDD)
#define ZGEMM  F77name(zgemm,ZGEMM) 
#define ZGEMV  F77name(zgemv,ZGEMV) 
