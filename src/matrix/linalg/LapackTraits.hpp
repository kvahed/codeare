#include "MLInterface.hpp"

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
        SPOTRF (&uplo, &n, a, &lda, &info);
    }
    
    inline static void 
	getrf (const int& m, const int& n, Type *a, const int& lda, int* ipiv, int& info) {
		SGETRF (&m, &n, a, &lda, ipiv, &info);
	}
    
    inline static void 
    potri (const char* uplo, int*n, void *a, int* lda, int*info) {
        SPOTRI (uplo, n, a, lda, info);
    }
    
    inline static void
    getri (const int& n, Type *a, const int& lda, int *ipiv, Type *work, const int& lwork,
           int& info) {
        SGETRI (&n, a, &lda, ipiv, work, &lwork, &info);
    }
    
    inline static void 
    gemv (const char& trans, const int& m, const int& n, const Type& alpha, const Type *a,
          const int& lda, const Type *x, const int& incx, const Type& beta, Type *y,
          const int& incy) {
        SGEMV (&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
    }
    
    inline static void 
    gemm (const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const Type& alpha, const Type *a, const int& lda, const Type *b, const int& ldb,
          const Type& beta, Type *c, const int& ldc) {
        SGEMM (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    }
    
    inline static double
    nrm2 (const int N, const Type *X, const int incX) {
        return SNRM2 (&N, X, &incX);
    }
    
    inline static void
    dot  (const int N, const Type *X, const int incX, const Type *Y, 
          const int incY, Type* res) {
        *res = SDOT (&N, X, &incX, Y, &incY);
    }
    
    inline static void
    dotc (const int N, const Type *X, const int incX, const Type *Y, 
          const int incY, Type* res) {
        *res = SDOT (&N, X, &incX, Y, &incY);
    }
    
    inline static void
    geev (const char& jvl, const char &jvr, const int& n, Type *a, const int& lda,
          CType *w, Type *vl, const int& ldvl, Type *vr, const int& ldvr, Type *work,
          const int& lwork, Type* rwork, int& info) {
        
        Type* dwr = (Type*) malloc (n * sizeof(Type));
        Type* dwi = (Type*) malloc (n * sizeof(Type));
        
        SGEEV (&jvl, &jvr, &n, a, &lda, dwr, dwi, vl, &ldvl, vr, &ldvr, work, &lwork,
               &info);
        for (int i = 0; i < n; i++)
            w[i] = CType (dwr[i], dwi[i]);
        
        free (dwr);
        free (dwi);
        
    }
    
    inline static void
    gels (const char& trans, const int& m, const int& n, const int& nrhs,
          Matrix<Type>& a, const int& lda, Matrix<Type>& b, const int& ldb,
          container<Type>& work, const int& lwork, int& info) {
        
        SGELS (&trans, &m, &n, &nrhs, &a[0], &lda, &b[0], &ldb, &work[0], &lwork, &info);
        
    }
    
    inline static void 
    gesdd (const char& jobz, const int& m, const int& n, Type *a, const int& lda, RType *s,
           Type* u, const int& ldu, Type *vt, const int& ldvt, Type* work, const int& lwork,
           RType *rwork, int* iwork, int& info) {
        SGESDD (&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);
    }
    
};


template <>
struct LapackTraits<double> {
    
    typedef double Type;
    typedef double RType;
    typedef std::complex<double> CType;
    
    inline static void 
    potrf (const char& uplo, const int& n, Type* a, const int& lda, int& info) {
		DPOTRF (&uplo, &n, a, &lda, &info);
	}
    
    inline static void 
	getrf (const int& m, const int& n, Type *a, const int& lda, int* ipiv, int& info) {
		DGETRF (&m, &n, a, &lda, ipiv, &info);
	}
    
    inline static void 
    potri (const char* uplo, int*n, void *a, int* lda, int*info) {
        DPOTRI (uplo, n, a, lda, info);
    }
    
    inline static void
	getri (const int& n, Type *a, const int& lda, int *ipiv, Type *work, const int& lwork,
		   int& info) {
		DGETRI (&n, a, &lda, ipiv, work, &lwork, &info);
	}
    
    inline static void 
    gemv (const char& trans, const int& m, const int& n, const Type& alpha, const Type *a,
          const int& lda, const Type *x, const int& incx, const Type& beta, Type *y,
          const int& incy) {
        DGEMV (&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
    }
    
    inline static void 
    gemm (const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const Type& alpha, const Type *a, const int& lda, const Type *b, const int& ldb,
          const Type& beta, Type *c, const int& ldc) {
        DGEMM (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    }

    inline static double
    nrm2 (const int N, const Type *X, const int incX) {
        return DNRM2 (&N, X, &incX);
    }
    
    inline static void
    dot  (const int N, const Type *X, const int incX, const Type *Y, 
          const int incY, Type* res) {
        *res = DDOT (&N, X, &incX, Y, &incY);
    }
    
    inline static void
    dotc (const int N, const Type *X, const int incX, const Type *Y, 
          const int incY, Type* res) {
        *res = DDOT (&N, X, &incX, Y, &incY);
    }

    inline static void
    geev (const char& jvl, const char &jvr, const int& n, Type *a, const int& lda,
          CType *w, Type *vl, const int& ldvl, Type *vr, const int& ldvr, Type *work,
          const int& lwork, Type* rwork, int& info) {
        
        Type* dwr = (Type*) malloc (n * sizeof(Type));
        Type* dwi = (Type*) malloc (n * sizeof(Type));
        
        DGEEV (&jvl, &jvr, &n, a, &lda, dwr, dwi, vl, &ldvl, vr, &ldvr, work, &lwork,
               &info);
        
        for (int i = 0; i < n; i++) {
            w[i] = CType (dwr[i], dwi[i]);
        }
        
        free (dwr);
        free (dwi);
        
    }
    
    inline static void
    gels (const char& trans, const int& m, const int& n, const int& nrhs,
          Matrix<Type>& a, const int& lda, Matrix<Type>& b, const int& ldb,
          container<Type>& work,
          const int& lwork, int& info) {
        
        DGELS (&trans, &m, &n, &nrhs, &a[0], &lda, &b[0], &ldb, &work[0], &lwork, &info);
        
    }

    inline static void 
    gesdd (const char& jobz, const int& m, const int& n, Type *a, const int& lda, RType *s,
           Type* u, const int& ldu, Type *vt, const int& ldvt, Type* work, const int& lwork,
           RType *rwork, int* iwork, int& info) {
        DGESDD (&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);
    }
    
};



template <>
struct LapackTraits<cxfl> {
    
    typedef cxfl Type;
    typedef float RType;
    typedef cxfl CType;
    
    inline static void 
    potrf (const char& uplo, const int& n, Type* a, const int& lda, int& info) {
		CPOTRF (&uplo, &n, a, &lda, &info);
	}
    
    inline static void 
	getrf (const int& m, const int& n, Type *a, const int& lda, int* ipiv, int& info) {
		CGETRF (&m, &n, a, &lda, ipiv, &info);
	}
    
    inline static void 
    potri (const char* uplo, int*n, void *a, int* lda, int*info) {
        CPOTRI (uplo, n, a, lda, info);
    }
    
    inline static void
	getri (const int& n, Type *a, const int& lda, int *ipiv, Type *work, const int& lwork,
		   int& info) {
		CGETRI (&n, a, &lda, ipiv, work, &lwork, &info);
	}
    
    inline static void 
    gemv (const char& trans, const int& m, const int& n, const Type& alpha, const Type *a,
          const int& lda, const Type *x, const int& incx, const Type& beta, Type *y,
          const int& incy) {
        CGEMV (&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
    }
    
    inline static void 
    gemm (const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const Type& alpha, const Type *a, const int& lda, const Type *b, const int& ldb,
          const Type& beta, Type *c, const int& ldc) {
        CGEMM (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    }
    
    inline static double
    nrm2 (const int N, const Type *X, const int incX) {
        return SCNRM2 (&N, X, &incX);
    }
    
    inline static void 
    dot  (const int N, const cxfl *X, const int incX, const cxfl *Y, 
          const int incY, cxfl* res) {
#if defined __APPLE__
        CDOTU (res, &N, X, &incX, Y, &incY);
#elif defined _MSC_VER
        *res = *X++ * *Y++;
        for (int i = 1; i < N; ++i)
            *res += *X++ * *Y++;
#else
        *res = CDOTU (&N, X, &incX, Y, &incY);
#endif
    }
    
    inline static void 
    dotc (const int N, const cxfl *X, const int incX, const cxfl *Y, 
          const int incY, cxfl* res) {
#if defined __APPLE__
        CDOTC (res, &N, X, &incX, Y, &incY);
#elif defined _MSC_VER
        *res = conj(*X++) * *Y++;
        for (int i = 1; i < N; ++i)
            *res += conj(*X++) * *Y++;
#else
        *res = CDOTC (&N, X, &incX, Y, &incY);
#endif
    }
    
    inline static void 
    geev (const char& jvl, const char& jvr, const int& n, Type *a, const int& lda,
          Type *w, Type *vl, const int& ldvl, Type *vr, const int& ldvr, Type *work,
          const int& lwork, RType* rwork, int& info) {
        CGEEV  (&jvl, &jvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork,
                rwork, &info);
    }
    
    inline static void
    gels (const char& trans, const int& m, const int& n, const int& nrhs,
          Matrix<Type>& a, const int& lda, Matrix<Type>& b, const int& ldb,
          container<Type>& work, const int& lwork, int& info) {
        
        CGELS (&trans, &m, &n, &nrhs, &a[0], &lda, &b[0], &ldb, &work[0], &lwork, &info);
        
    }
    
    inline static void 
    gesdd (const char& jobz, const int& m, const int& n, Type *a, const int& lda, RType *s,
           Type* u, const int& ldu, Type *vt, const int& ldvt, Type* work, const int& lwork,
           RType *rwork, int* iwork, int& info) {
        CGESDD (&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info);
    }
    
};



template <>
struct LapackTraits<cxdb> {
    
    typedef cxdb Type;
    typedef double RType;
    typedef cxdb CType;
    
    inline static void 
    potrf (const char& uplo, const int& n, Type* a, const int& lda, int& info) {
		ZPOTRF (&uplo, &n, a, &lda, &info);
	}
    
    inline static void 
    getrf (const int& m, const int& n, Type *a, const int& lda, int* ipiv, int& info) {
        ZGETRF (&m, &n, a, &lda, ipiv, &info);
    }
    
    inline static void 
    potri (const char* uplo, int*n, void *a, int* lda, int*info) {
        ZPOTRI (uplo, n, a, lda, info);
    }
    
    inline static void
	getri (const int& n, Type *a, const int& lda, int *ipiv, Type *work, const int& lwork,
		   int& info) {
		ZGETRI (&n, a, &lda, ipiv, work, &lwork, &info);
	}
    
    inline static void 
    gemv (const char& trans, const int& m, const int& n, const Type& alpha, const Type *a,
          const int& lda, const Type *x, const int& incx, const Type& beta, Type *y,
          const int& incy) {
        ZGEMV (&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
    }
    
    inline static void 
    gemm (const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const Type& alpha, const Type *a, const int& lda, const Type *b, const int& ldb,
          const Type& beta, Type *c, const int& ldc) {
        ZGEMM (&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    }
    
    inline static double
    nrm2 (const int N, const Type *X, const int incX) {
        return DZNRM2 (&N, X, &incX);
    }
    
    inline static void 
    dot  (const int N, const cxdb *X, const int incX, const cxdb *Y, 
          const int incY, cxdb* res) {
#if defined __APPLE__
        ZDOTU (res, &N, X, &incX, Y, &incY);
#elif defined _MSC_VER
        *res = *X++ * *Y++;
        for (int i = 1; i < N; ++i)
            *res += *X++ * *Y++;
#else
        *res = ZDOTU (&N, X, &incX, Y, &incY);
#endif
    }
    
    inline static void 
    dotc (const int N, const cxdb *X, const int incX, const cxdb *Y, 
          const int incY, cxdb* res) {
#if defined __APPLE__
        ZDOTC (res, &N, X, &incX, Y, &incY);
#elif defined _MSC_VER
        *res = conj(*X++) * *Y++;
        for (int i = 1; i < N; ++i)
            *res += conj(*X++) * *Y++;
#else
        *res = ZDOTC (&N, X, &incX, Y, &incY);
#endif
    }
    
    inline static void 
    geev (const char& jvl, const char& jvr, const int& n, Type *a, const int& lda,
          Type *w, Type *vl, const int& ldvl, Type *vr, const int& ldvr, cxdb *work,
          const int& lwork, double* rwork, int& info) {
        ZGEEV  (&jvl, &jvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork,
                rwork, &info);
    }
    
    inline static void
    gels (const char& trans, const int& m, const int& n, const int& nrhs,
          Matrix<Type>& a, const int& lda, Matrix<Type>& b, const int& ldb,
          container<Type>& work, const int& lwork, int& info) {
        
        ZGELS (&trans, &m, &n, &nrhs, &a[0], &lda, &b[0], &ldb, &work[0], &lwork, &info);
        
    }
    
    inline static void 
    gesdd (const char& jobz, const int& m, const int& n, Type *a, const int& lda, RType *s,
           Type* u, const int& ldu, Type *vt, const int& ldvt, Type* work, const int& lwork,
           RType *rwork, int* iwork, int& info) {
        ZGESDD (&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info);
    }
    
    
};


template <>
struct LapackTraits<unsigned char> {

    typedef unsigned char Type;
    typedef unsigned char RType;
    typedef unsigned char CType;

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
          Matrix<Type>& a, const int& lda, Matrix<Type>& b, const int& ldb,
          container<Type>& work, const int& lwork, int& info) {}

    inline static void
    gesdd (const char& jobz, const int& m, const int& n, Type *a, const int& lda, RType *s,
           Type* u, const int& ldu, Type *vt, const int& ldvt, Type* work, const int& lwork,
           RType *rwork, int* iwork, int& info) {}

};

template <>
struct LapackTraits<size_t> {

    typedef size_t Type;
    typedef size_t RType;
    typedef size_t CType;

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
          Matrix<Type>& a, const int& lda, Matrix<Type>&, const int& ldb,
          container<Type>& work, const int& lwork, int& info) {}

    inline static void
    gesdd (const char& jobz, const int& m, const int& n, Type *a, const int& lda, RType *s,
           Type* u, const int& ldu, Type *vt, const int& ldvt, Type* work, const int& lwork,
           RType *rwork, int* iwork, int& info) {}

};
