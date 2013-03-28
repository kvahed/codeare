#include "Complex.hpp"

#ifndef __SCALAPACK_TRAITS_HPP__
#define __SCALAPACK_TRAITS_HPP__

extern "C" {
	
	// BLACS
	void Cblacs_pinfo (int* mypnum, int* nprocs);
	void Cblacs_get (int context, int request, int* value);
	int  Cblacs_gridinit (int* context, char* order, int np_row, int np_col);
	void Cblacs_gridinfo (int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
	void Cblacs_gridexit (int context);
	void Cblacs_exit (int error_code);
	void Cblacs_barrier (int context, char* scope);
	
    void Cdgerv2d (int, int, int, double*, int, int, int);
    void Cdgesd2d (int, int, int, double*, int, int, int);

	// Initialise descriptor vector 
    void descinit_  (int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, 
					 int *icsrc, int *ictxt, int *lld, int *info); 
	
	// GLobal index from local
	int  indxl2g_   (int *lidx , int *nb, int *iproc, int *isrcproc, int *nprocs);
	
	// Local grid 
    int  numroc_    (int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);

	// SVD
    void psgesvd_   (const char* jobu, const char* jobvt, const int* m, const int* n,
                     float* A, const int* ia, const int* ja, const int* desca,
                     float* S, float* U, const int* iu, const int* ju,
                     const int* descu, float* VT, const int* ivt, const int* jvt,
                     const int* descvt, float* work, const int* lwork, int* info);
    void pdgesvd_   (const char* jobu, const char* jobvt, const int* m, const int* n,
                     double* A, const int* ia, const int* ja, const int* desca,
                     double* S, double* U, const int* iu, const int* ju,
                     const int* descu, double* VT, const int* ivt, const int* jvt,
                     const int* descvt, double* work, const int* lwork, int* info);
    void pcgesvd_   (const char* jobu, const char* jobvt, const int* m, const int* n,
                     cxfl* A, const int* ia, const int* ja, const int* desca,
                     float* S, cxfl* U, const int* iu, const int* ju,
                     const int* descu, cxfl* VT, const int* ivt, const int* jvt,
                     const int* descvt, cxfl* work, const int* lwork, float* rwork,
                     int* info);
    void pzgesvd_   (const char* jobu, const char* jobvt, const int* m, const int* n,
                     cxdb* A, const int* ia, const int* ja, const int* desca,
                     double* S, cxdb* U, const int* iu, const int* ju,
                     const int* descu, cxdb* VT, const int* ivt, const int* jvt,
                     const int* descvt, cxdb* work, const int* lwork, double* rwork,
                     int* info);
	
	// Print
    void pzlaprnt_  (int *m, int* n, const cxdb* A, int* ia, int* ja, 
					 const int* descA, int* irprnt, int* icprnt, 
					 const char* cmatnm, int* nout, cxdb* WORK, int len);
    void pclaprnt_  (int *m, int* n, const cxfl* A, int* ia, int* ja, 
					 const int* descA, int* irprnt, int* icprnt, 
					 const char* cmatnm, int* nout, cxfl* WORK, int len);
    void pdlaprnt_  (int *m, int* n, const double* A, int* ia, int* ja, 
					 const int* descA, int* irprnt, int* icprnt, 
					 const char* cmatnm, int* nout, double* WORK, int len);
    void pslaprnt_  (int *m, int* n, const float* A, int* ia, int* ja, 
					 const int* descA, int* irprnt, int* icprnt, 
					 const char* cmatnm, int* nout, float* WORK, int len);
	
	// File IO
	void pzlawrite_ (char* fname, int* m, int* n, cxdb* A, int* ia, int* ja, 
					 int* descA,	int* irwr, int* icwr, cxdb* WORK);
	void pclawrite_ (char* fname, int* m, int* n, cxdb* A, int* ia, int* ja, 
					 int* descA,	int* irwr, int* icwr, cxdb* WORK);
	
	// Solve LSQR under / over determined ||Ax-b||
    void psgels_    (const char* trans, const int* m, const int* n, const int* nrhs,
                     float* A, const int* ia, const int* ja, const int* desca,
                     float* B, const int* ib, const int* jb, const int* descb,
                     float* work, const int* lwork, int* info);
    void pdgels_    (const char* trans, const int* m, const int* n, const int* nrhs,
                     double* A, const int* ia, const int* ja, const int* desca,
                     double* B, const int* ib, const int* jb, const int* descb,
                     double* work, const int* lwork, int* info);
    void pcgels_    (const char* trans, const int* m, const int* n, const int* nrhs,
                     cxfl* A, const int* ia, const int* ja, const int* desca,
                     cxfl* B, const int* ib, const int* jb, const int* descb,
                     cxfl* work, const int* lwork, int* info);
    void pzgels_    (const char* trans, const int* m, const int* n, const int* nrhs,
                     cxdb* A, const int* ia, const int* ja, const int* desca,
                     cxdb* B, const int* ib, const int* jb, const int* descb,
                     cxdb* work, const int* lwork, int* info);

	// MV multiplication (y <- alpha * sub(op(A)) * sub(x) + beta * sub(y))
	void psgemv_    (const char* trans, const int* m, const int* n, const float* alpha,
                     const float* A, const int* ia, const int* ja, const int* desca,
                     const float* X, const int* ix, const int* jx, const int* descx,
                     const int* incx, const float* beta, float* y, const int* iy,
                     const int* jy, const int* descY, const int* incy);
	void pdgemv_    (const char* trans, const int* m, const int* n, const double* alpha,
                     const double* A, const int* ia, const int* ja, const int* desca,
                     const double* X, const int* ix, const int* jx, const int* descx,
                     const int* incx, const double* beta, double* y, const int* iy,
                     const int* jy, const int* descY, const int* incy);
	void pcgemv_    (const char* trans, const int* m, const int* n, const cxfl* alpha,
                     const cxfl* A, const int* ia, const int* ja, const int* desca,
                     const cxfl* X, const int* ix, const int* jx, const int* descx,
                     const int* incx, const cxfl* beta, cxfl* y, const int* iy,
                     const int* jy, const int* descY, const int* incy);
	void pzgemv_    (const char* trans, const int* m, const int* n, const cxdb* alpha,
                     const cxdb* A, const int* ia, const int* ja, const int* desca,
                     const cxdb* X, const int* ix, const int* jx, const int* descx,
                     const int* incx, const cxdb* beta, cxdb* y, const int* iy,
                     const int* jy, const int* descY, const int* incy);
	
	// MM multiplication (C <- alpha * sub(op(A)) * sub(op(B))) + beta * C
	void psgemm_    (const char* transa, const char* transb, const int *m,
                     const int *n, const int *k, float *alpha, const float *A,
                     const int *ia, const int *ja, int * desca, const float *B,
                     const int *ib, const int *jb, const int * descb,
                     const float *beta, float *C, const int* ic, const int* jc,
                     const int *descc);
	void pdgemm_    (const char* transa, const char* transb, const int *m,
                     const int *n, const int *k, double *alpha, const double *A,
                     const int *ia, const int *ja, int * desca, const double *B,
                     const int *ib, const int *jb, const int * descb,
                     const double *beta, double *C, const int* ic, const int* jc,
                     const int *descc);
	void pcgemm_    (const char* transa, const char* transb, const int *m,
                     const int *n, const int *k, cxfl *alpha, const cxfl *A,
                     const int *ia, const int *ja, int * desca, const cxfl *B,
                     const int *ib, const int *jb, const int * descb,
                     const cxfl *beta, cxfl *C, const int* ic, const int* jc,
                     const int *descc);
	void pzgemm_    (const char* transa, const char* transb, const int *m,
                     const int *n, const int *k, cxdb *alpha, const cxdb *A,
                     const int *ia, const int *ja, int * desca, const cxdb *B,
                     const int *ib, const int *jb, const int * descb,
                     const cxdb *beta, cxdb *C, const int* ic, const int* jc,
                     const int *descc);


	// Cholesky factorisation of an NxN [hermitian] positive definite matrix
	void pcpotrf_ (const char* uplo, const int& n, cxfl* A, const int& ia, 
				   const int& ja, const int* descA, int& info);
	void pzpotrf_ (const char* uplo, const int& n, cxdb* A, const int& ia, 
				   const int& ja, const int* descA, int& info);
	void pspotrf_ (const char* uplo, const int& n, float* A, const int& ia, 
				   const int& ja, const int* descA, int& info);
	void pdpotrf_ (const char* uplo, const int& n, double* A, const int& ia, 
				   const int& ja, const int* descA, int& info);

	// LU factorization of a general MxN matrix
	void psgetrf_ (const int& m, const int& n, float* a, const int* lda, 
				   int* ipiv, const int& info);
	void pdgetrf_ (const int& m, const int& n, double* a, const int* lda, 
				   int* ipiv, const int& info);
	void pcgetrf_ (const int& m, const int& n, cxfl* a, const int* lda, 
				   int* ipiv, const int& info);
	void pzgetrf_ (const int& m, const int& n, cxdb* a, const int* lda, 
				   int* ipiv, const int& info);
	
	// Inverse of a complex Hermitian pos def square mat with cpotrf/cpptrf
	void pcpotri_ (const char* uplo, const int& n, cxfl *a, const int* lda, 
				   int& info);
	void pdpotri_ (const char* uplo, const int& n, void *a, const int* lda, 
				   int& info);
	void pzpotri_ (const char* uplo, const int& n, void *a, const int* lda, 
				   int& info);
	void pspotri_ (const char* uplo, const int& n, void *a, const int* lda, 
				   int& info);

}


// C++ convenience
inline int indxl2g (int idl, int nb, int iproc, int isrcproc, int nprocs) {
	int fortidl = idl + 1;
	return indxl2g_ (&fortidl, &nb, &iproc, &isrcproc, &nprocs) - 1;
}




template<class T>
struct ScalapackTraits {};



template<>
struct ScalapackTraits<cxfl> {
	
	typedef cxfl Type;
	
	inline static void
	gemv  (char* trans, int* m, int* n, cxfl* alpha, cxfl* A, int* ia, int* ja, 
		   int* descA, cxfl* X, int* ix, int* jx, int* descX, int* incx, 
		   cxfl* beta, cxfl* Y, int* iy, int* jy, int* descY, int* incy) {
		pcgemv_  (trans, m, n, alpha, A, ia, ja, descA, X, ix, jx, descX, incx, 
				  beta, Y, iy, jy, descY, incy);
	}
	
	inline static void 
	gemm  (char* transA, char* transB, int *m, int *n, int *k,	cxfl *alpha, 
		   cxfl *A, int *ia, int *ja, int *descA, cxfl *B, int *ib, int *jb, 
		   int *descB, cxfl *beta, cxfl *C, int * ic, int * jc, int *descC) {
		pcgemm_  (transA, transB, m, n, k, alpha, A, ia, ja, descA, B, ib, jb, 
				  descB, beta, C, ic, jc, descC);
	}

    inline static void 
	gesvd (char* jbu, char* jbvt, int* m, int* n, cxfl* A, int* ia, int* ja, 
		   int* descA, float* s, cxfl* U, int* iu, int* ju, int* descU, 
		   cxfl* VT, int* ivt, int* jvt, int* descVT, cxfl* WORK, int* lwork, 
		   float* rwork, int* info) {
		pcgesvd_ (jbu, jbvt, m, n, A, ia, ja, descA, s, U, iu, ju, descU, VT, ivt,
				  jvt, descVT, WORK, lwork, rwork, info);
	}

    inline static void 
	gels  (char* trans, int* m, int* n, int* nrhs, cxfl* A, int* ia, int* ja, 
		   int* descA, cxfl* B, int* ib, int* jb, int* descB, cxfl* WORK, 
		   int* lwork, int* info) {
		pcgels_  (trans, m, n, nrhs, A, ia, ja, descA, B, ib, jb, descB, WORK,
				  lwork, info);
	}


};

template<>
struct ScalapackTraits<cxdb> {
	
	typedef cxdb Type;
	
	inline static void
	gemv (char* trans, int* m, int* n, cxdb* alpha, cxdb* A, int* ia, int* ja, 
		  int* descA, cxdb* X, int* ix, int* jx, int* descX, int* incx, 
		  cxdb* beta, cxdb* Y, int* iy, int* jy, int* descY, int* incy) {
		pzgemv_ (trans, m, n, alpha, A, ia, ja, descA, X, ix, jx, descX, incx, 
				 beta, Y, iy, jy, descY, incy);
	}
	
	inline static void 
	gemm  (char* transA, char* transB, int *m, int *n, int *k,	cxdb *alpha, 
		   cxdb *A, int *ia, int *ja, int *descA, cxdb *B, int *ib, int *jb, 
		   int *descB, cxdb *beta, cxdb *C, int * ic, int * jc, int *descC) {
		pzgemm_  (transA, transB, m, n, k, alpha, A, ia, ja, descA, B, ib, jb, 
				  descB, beta, C, ic, jc, descC);
	}

    inline static void 
	gesvd (char* jbu, char* jbvt, int* m, int* n, cxdb* A, int* ia, int* ja, 
		   int* descA, double* s, cxdb* U, int* iu, int* ju, int* descU, 
		   cxdb* VT, int* ivt, int* jvt, int* descVT, cxdb* WORK, int* lwork, 
		   double* rwork, int* info) {
		pzgesvd_ (jbu, jbvt, m, n, A, ia, ja, descA, s, U, iu, ju, descU, VT, ivt,
				  jvt, descVT, WORK, lwork, rwork, info);
	}	

    inline static void 
	gels  (char* trans, int* m, int* n, int* nrhs, cxdb* A, int* ia, int* ja, 
		   int* descA, cxdb* B, int* ib, int* jb, int* descB, cxdb* WORK, 
		   int* lwork, int* info) {
		pzgels_ (trans, m, n, nrhs, A, ia, ja, descA, B, ib, jb, descB, WORK,
				 lwork, info);
	}

};


template<>
struct ScalapackTraits<double> {

	typedef double Type;

	inline static void 
	pxlaprnt (int *m, int* n, const double* A, int* ia, int* ja, const int* descA, 
			  int* irprnt, int* icprnt, const char* cmatnm, int* nout, 
			  double* WORK, int len) {
		char scope = 'A';
		Cblacs_barrier (descA[1], &scope);
		pdlaprnt_ (m, n, A, ia, ja, descA, irprnt, icprnt, cmatnm, nout, WORK, len);
	}

};

template<>
struct ScalapackTraits<float> {

	typedef float Type;

	inline static void 
	pxlaprnt (int *m, int* n, const float* A, int* ia, int* ja, const int* descA, 
			  int* irprnt, int* icprnt, const char* cmatnm, int* nout, 
			  float* WORK, int len) {
		char scope = 'A';
		Cblacs_barrier (descA[1], &scope);
		pslaprnt_ (m, n, A, ia, ja, descA, irprnt, icprnt, cmatnm, nout, WORK, len);
	}

};

#endif
