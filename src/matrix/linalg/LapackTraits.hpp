#include "MLInterface.hpp"

#include <cblas.h>
#include <string.h>

template <class T>
struct LapackTraits {};


template <>
struct LapackTraits<float> {
    
    typedef float Type;
    typedef float RType;
    typedef std::complex<float> CType;
    
    
    inline static void 
    potrf (const char& uplo, const int& n, Type* a, const int& lda, int& info) {
#ifdef USE_ACML
#else
        SPOTRF (&uplo, &n, a, &lda, &info);
#endif
    }
    
    inline static void 
	getrf (const int& m, const int& n, Type *a, const int& lda, int* ipiv, int& info) {
#ifdef USE_ACML
#else
		SGETRF (&m, &n, a, &lda, ipiv, &info);
#endif
	}
    
    inline static void 
    potri (const char* uplo, int*n, void *a, int* lda, int*info) {
#ifdef USE_ACML
#else
        SPOTRI (uplo, n, a, lda, info);
#endif
    }
    
    inline static void
    getri (const int& n, Type *a, const int& lda, int *ipiv, Type *work, const int& lwork,
           int& info) {
#ifdef USE_ACML
#else
        SGETRI (&n, a, &lda, ipiv, work, &lwork, &info);
#endif
    }
    
    inline static void 
    gemv (const char& trans, const int& m, const int& n, const Type& alpha, const Type *a,
          const int& lda, const Type *x, const int& incx, const Type& beta, Type *y,
          const int& incy) {
#ifdef USE_ACML
#else
        SGEMV (&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endif
    }
    
    inline static void 
    gemm (const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const Type& alpha, const Type *a, const int& lda, const Type *b, const int& ldb,
          const Type& beta, Type *c, const int& ldc) {
#ifdef USE_ACML
#else
        SGEMM (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
    }
    
    inline static Type 
    nrm2 (const int N, const Type *X, const int incX) {
        return cblas_snrm2 (N, X, incX);
    }
    
    inline static void
    dot  (const int N, const Type *X, const int incX, const Type *Y, 
          const int incY, Type* res) {
        *res = cblas_sdot (N, X, incX, Y, incY);
    }
    
    inline static void
    dotc (const int N, const Type *X, const int incX, const Type *Y, 
          const int incY, Type* res) {
        *res = cblas_sdot (N, X, incX, Y, incY);
    }
    
    inline static void
    geev (const char& jvl, const char &jvr, const int& n, Type *a, const int& lda,
          CType *w, Type *vl, const int& ldvl, Type *vr, const int& ldvr, Type *work,
          const int& lwork, Type* rwork, int& info) {
        
        Type* dwr = (Type*) malloc (n * sizeof(Type));
        Type* dwi = (Type*) malloc (n * sizeof(Type));
        
#ifdef USE_ACML
#else
        SGEEV (&jvl, &jvr, &n, a, &lda, dwr, dwi, vl, &ldvl, vr, &ldvr, work, &lwork,
               &info);
#endif
        
        for (size_t i = 0; i < n; i++) {
            w[i] = CType (dwr[i], dwi[i]);
        }
        
        free (dwr);
        free (dwi);
        
    }
    
    inline static void
    gels (const char& trans, const int& m, const int& n, const int& nrhs,
          Type* a, const int& lda, Type* b, const int& ldb, Type* work,
          const int& lwork, int& info) {
        
        SGELS (&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
        
    }
    
    inline static void 
    gesdd (const char& jobz, const int& m, const int& n, Type *a, const int& lda, RType *s,
           Type* u, const int& ldu, Type *vt, const int& ldvt, Type* work, const int& lwork,
           RType *rwork, int* iwork, int& info) {
#ifdef USE_ACML
#else
        SGESDD (&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);
#endif
    }
    
};


template <>
struct LapackTraits<double> {
    
    typedef double Type;
    typedef double RType;
    typedef std::complex<double> CType;
    
    inline static void 
    potrf (const char& uplo, const int& n, Type* a, const int& lda, int& info) {
#ifdef USE_ACML
#else
		DPOTRF (&uplo, &n, a, &lda, &info);
#endif
	}
    
    inline static void 
	getrf (const int& m, const int& n, Type *a, const int& lda, int* ipiv, int& info) {
#ifdef USE_ACML
#else
		DGETRF (&m, &n, a, &lda, ipiv, &info);
#endif
	}
    
    inline static void 
    potri (const char* uplo, int*n, void *a, int* lda, int*info) {
#ifdef USE_ACML
#else
        DPOTRI (uplo, n, a, lda, info);
#endif
    }
    
    inline static void
	getri (const int& n, Type *a, const int& lda, int *ipiv, Type *work, const int& lwork,
		   int& info) {
#ifdef USE_ACML
#else
		DGETRI (&n, a, &lda, ipiv, work, &lwork, &info);
#endif
	}
    
    inline static void 
    gemv (const char& trans, const int& m, const int& n, const Type& alpha, const Type *a,
          const int& lda, const Type *x, const int& incx, const Type& beta, Type *y,
          const int& incy) {
#ifdef USE_ACML
#else
        DGEMV (&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endif
    }
    
    inline static void 
    gemm (const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const Type& alpha, const Type *a, const int& lda, const Type *b, const int& ldb,
          const Type& beta, Type *c, const int& ldc) {
#ifdef USE_ACML
#else
        DGEMM (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
    }

    inline static Type 
    nrm2 (const int N, const Type *X, const int incX) {
        return cblas_dnrm2 (N, X, incX);
    }
    
    inline static void
    dot  (const int N, const Type *X, const int incX, const Type *Y, 
          const int incY, Type* res) {
        *res = cblas_ddot (N, X, incX, Y, incY);
    }
    
    inline static void
    dotc (const int N, const Type *X, const int incX, const Type *Y, 
          const int incY, Type* res) {
        *res = cblas_ddot (N, X, incX, Y, incY);
    }

    inline static void
    geev (const char& jvl, const char &jvr, const int& n, Type *a, const int& lda,
          CType *w, Type *vl, const int& ldvl, Type *vr, const int& ldvr, Type *work,
          const int& lwork, Type* rwork, int& info) {
        
        Type* dwr = (Type*) malloc (n * sizeof(Type));
        Type* dwi = (Type*) malloc (n * sizeof(Type));
        
#ifdef USE_ACML
#else
        DGEEV (&jvl, &jvr, &n, a, &lda, dwr, dwi, vl, &ldvl, vr, &ldvr, work, &lwork,
               &info);
#endif
        
        for (size_t i = 0; i < n; i++) {
            w[i] = CType (dwr[i], dwi[i]);
        }
        
        free (dwr);
        free (dwi);
        
    }
    
    inline static void
    gels (const char& trans, const int& m, const int& n, const int& nrhs,
          Type* a, const int& lda, Type* b, const int& ldb, Type* work,
          const int& lwork, int& info) {
        
        DGELS (&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
        
    }

    inline static void 
    gesdd (const char& jobz, const int& m, const int& n, Type *a, const int& lda, RType *s,
           Type* u, const int& ldu, Type *vt, const int& ldvt, Type* work, const int& lwork,
           RType *rwork, int* iwork, int& info) {
#ifdef USE_ACML
#else
        DGESDD (&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);
#endif
    }
    
};



template <>
struct LapackTraits<cxfl> {
    
    typedef cxfl Type;
    typedef float RType;
    typedef cxfl CType;
    
    inline static void 
    potrf (const char& uplo, const int& n, Type* a, const int& lda, int& info) {
#ifdef USE_ACML
#else
		CPOTRF (&uplo, &n, a, &lda, &info);
#endif
	}
    
    inline static void 
	getrf (const int& m, const int& n, Type *a, const int& lda, int* ipiv, int& info) {
#ifdef USE_ACML
#else
		CGETRF (&m, &n, a, &lda, ipiv, &info);
#endif
	}
    
    inline static void 
    potri (const char* uplo, int*n, void *a, int* lda, int*info) {
#ifdef USE_ACML
#else
        CPOTRI (uplo, n, a, lda, info);
#endif
    }
    
    inline static void
	getri (const int& n, Type *a, const int& lda, int *ipiv, Type *work, const int& lwork,
		   int& info) {
#ifdef USE_ACML
#else
		CGETRI (&n, a, &lda, ipiv, work, &lwork, &info);
#endif
	}
    
    inline static void 
    gemv (const char& trans, const int& m, const int& n, const Type& alpha, const Type *a,
          const int& lda, const Type *x, const int& incx, const Type& beta, Type *y,
          const int& incy) {
#ifdef USE_ACML
#else
        CGEMV (&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endif
    }
    
    inline static void 
    gemm (const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const Type& alpha, const Type *a, const int& lda, const Type *b, const int& ldb,
          const Type& beta, Type *c, const int& ldc) {
#ifdef USE_ACML
#else
        CGEMM (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
    }
    
    inline static Type
    nrm2 (const int N, const Type *X, const int incX) {
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
    geev (const char& jvl, const char& jvr, const int& n, Type *a, const int& lda,
          Type *w, Type *vl, const int& ldvl, Type *vr, const int& ldvr, Type *work,
          const int& lwork, RType* rwork, int& info) {
#ifdef USE_ACML
#else
        CGEEV  (&jvl, &jvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork,
                rwork, &info);
#endif
    }
    
    inline static void
    gels (const char& trans, const int& m, const int& n, const int& nrhs,
          Type* a, const int& lda, Type* b, const int& ldb, Type* work,
          const int& lwork, int& info) {
        
        CGELS (&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
        
    }
    
    inline static void 
    gesdd (const char& jobz, const int& m, const int& n, Type *a, const int& lda, RType *s,
           Type* u, const int& ldu, Type *vt, const int& ldvt, Type* work, const int& lwork,
           RType *rwork, int* iwork, int& info) {
#ifdef USE_ACML
#else
        CGESDD (&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info);
#endif
    }
    
};



template <>
struct LapackTraits<cxdb> {
    
    typedef cxdb Type;
    typedef double RType;
    typedef cxdb CType;
    
    inline static void 
    potrf (const char& uplo, const int& n, Type* a, const int& lda, int& info) {
#ifdef USE_ACML
#else
		ZPOTRF (&uplo, &n, a, &lda, &info);
#endif
	}
    
    inline static void 
    getrf (const int& m, const int& n, Type *a, const int& lda, int* ipiv, int& info) {
#ifdef USE_ACML
#else
        ZGETRF (&m, &n, a, &lda, ipiv, &info);
#endif
    }
    
    inline static void 
    potri (const char* uplo, int*n, void *a, int* lda, int*info) {
#ifdef USE_ACML
#else
        ZPOTRI (uplo, n, a, lda, info);
#endif
    }
    
    inline static void
	getri (const int& n, Type *a, const int& lda, int *ipiv, Type *work, const int& lwork,
		   int& info) {
#ifdef USE_ACML
#else
		ZGETRI (&n, a, &lda, ipiv, work, &lwork, &info);
#endif
	}
    
    inline static void 
    gemv (const char& trans, const int& m, const int& n, const Type& alpha, const Type *a,
          const int& lda, const Type *x, const int& incx, const Type& beta, Type *y,
          const int& incy) {
#ifdef USE_ACML
#else
        ZGEMV (&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endif
    }
    
    inline static void 
    gemm (const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const Type& alpha, const Type *a, const int& lda, const Type *b, const int& ldb,
          const Type& beta, Type *c, const int& ldc) {
#ifdef USE_ACML
#else
        ZGEMM (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
    }
    
    inline static Type
    nrm2 (const int N, const Type *X, const int incX) {
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
    geev (const char& jvl, const char& jvr, const int& n, Type *a, const int& lda,
          Type *w, Type *vl, const int& ldvl, Type *vr, const int& ldvr, cxdb *work,
          const int& lwork, double* rwork, int& info) {
#ifdef USE_ACML
#else
        ZGEEV  (&jvl, &jvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork,
                rwork, &info);
#endif
    }
    
    inline static void
    gels (const char& trans, const int& m, const int& n, const int& nrhs,
          Type* a, const int& lda, Type* b, const int& ldb, Type* work,
          const int& lwork, int& info) {
        
        ZGELS (&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
        
    }
    
    inline static void 
    gesdd (const char& jobz, const int& m, const int& n, Type *a, const int& lda, RType *s,
           Type* u, const int& ldu, Type *vt, const int& ldvt, Type* work, const int& lwork,
           RType *rwork, int* iwork, int& info) {
#ifdef USE_ACML
#else
        ZGESDD (&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info);
#endif
    }
    
    
};


template <>
struct LapackTraits<bool> {

    typedef bool Type;
    typedef bool RType;
    typedef bool CType;


    inline static void
    potrf (const char& uplo, const int& n, Type* a, const int& lda, int& info) {}

    inline static void
	getrf (const int& m, const int& n, Type *a, const int& lda, int* ipiv, int& info) {}

    inline static void
    potri (const char* uplo, int*n, void *a, int* lda, int*info) {}

    inline static void
    getri (const int& n, Type *a, const int& lda, int *ipiv, Type *work, const int& lwork,
           int& info) {}

    inline static void
    gemv (const char& trans, const int& m, const int& n, const Type& alpha, const Type *a,
          const int& lda, const Type *x, const int& incx, const Type& beta, Type *y,
          const int& incy) {}

    inline static void
    gemm (const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const Type& alpha, const Type *a, const int& lda, const Type *b, const int& ldb,
          const Type& beta, Type *c, const int& ldc) {}

    inline static Type
    nrm2 (const int N, const Type *X, const int incX) {
    	return false;
    }

    inline static Type
    dot  (const int N, const Type *X, const int incX, const Type *Y,
          const int incY, void* res) {
    	return false;
    }

    inline static Type
    dotc (const int N, const Type *X, const int incX, const Type *Y,
          const int incY, void* res) {
    	return false;
    }

    inline static void
    geev (const char& jvl, const char &jvr, const int& n, Type *a, const int& lda,
          CType *w, Type *vl, const int& ldvl, Type *vr, const int& ldvr, Type *work,
          const int& lwork, Type* rwork, int& info) {}

    inline static void
    gels (const char& trans, const int& m, const int& n, const int& nrhs,
          Type* a, const int& lda, Type* b, const int& ldb, Type* work,
          const int& lwork, int& info) {}

    inline static void
    gesdd (const char& jobz, const int& m, const int& n, Type *a, const int& lda, RType *s,
           Type* u, const int& ldu, Type *vt, const int& ldvt, Type* work, const int& lwork,
           RType *rwork, int* iwork, int& info) {}

};
