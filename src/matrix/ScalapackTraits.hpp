#include "Complex.hpp"

#define  MPICH_IGNORE_CXX_SEEK
#include <mpi.h>

struct grid_dims {
	
	int np; /**< @brief # of processes */
	int rk; /**< @brief my rank        */
	int ct; /**< @brief context        */
	int nr; /**< @brief # of proc rows */
	int nc; /**< @brief # of proc columns */
	int mr; /**< @brief my row # */
	int mc; /**< @brief my col # */ 
	char order; /**< @brief row/col major */

};

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
    void pzgesvd_   (char* jbu, char* jbvt, int* m, int* n, cxdb* A, int* ia, 
					 int* ja, int* descA, double* s, cxdb* U, int* iu, int* ju, 
					 int* descU, cxdb* VT, int* ivt, int* jvt, int* descVT, 
					 cxdb* WORK, int* lwork, double* rwork, int* info);
    void pcgesvd_   (char* jbu, char* jbvt, int* m, int* n, cxfl* A, int* ia, 
					 int* ja, int* descA, double* s, cxfl* U, int* iu, int* ju, 
					 int* descU, cxfl* VT, int* ivt, int* jvt, int* descVT, 
					 cxfl* WORK, int* lwork, double* rwork, int* info);
	
	// Print
    void pzlaprnt_  (int *m, int* n, cxdb* A, int* ia, int* ja, int* descA, 
					 int* irprnt, int* icprnt, char* cmatnm, int* nout, 
					 cxdb* WORK, int len);
    void pclaprnt_  (int *m, int* n, cxfl* A, int* ia, int* ja, int* descA, 
					 int* irprnt, int* icprnt, char* cmatnm, int* nout, 
					 cxfl* WORK, int len);

	// File IO
	void pzlawrite_ (char* fname, int* m, int* n, cxdb* A, int* ia, int* ja, 
					 int* descA,	int* irwr, int* icwr, cxdb* WORK);
	void pclawrite_ (char* fname, int* m, int* n, cxdb* A, int* ia, int* ja, 
					 int* descA,	int* irwr, int* icwr, cxdb* WORK);
	
	// Solve LSQR under / over determined ||Ax-b||
    void pzgels_    (char* trans, int* m, int* n, int* nrhs, cxdb* A, int* ia, 
					 int* ja, int* descA, cxdb* B, int* ib, int* jb, int* descB, 
					 cxdb* WORK, int* lwork, int* info);
    void pcgels_    (char* trans, int* m, int* n, int* nrhs, cxfl* A, int* ia, 
					 int* ja, int* descA, cxfl* B, int* ib, int* jb, int* descB, 
					 cxfl* WORK, int* lwork, int* info);

	// MV multiplication (y <- alpha * sub(op(A)) * sub(x) + beta * sub(y))
	void pzgemv_    (char* trans, int* m, int* n, cxdb* alpha, cxdb* A, 
					 int* ia, int* ja, int* descA, cxdb* x, int* ix, int* jx, 
					 int* descX, int* incx, cxdb* beta, cxdb* y, int* iy, int* jy, 
					 int* descY, int* incy);
	void pcgemv_    (char* trans, int* m, int* n, cxfl* alpha, cxfl* A, 
					 int* ia, int* ja, int* descA, cxfl* x, int* ix, int* jx, 
					 int* descX, int* incx, cxfl* beta, cxfl* y, int* iy, int* jy, 
					 int* descY, int* incy);
	
	// MM multiplication (C <- alpha * sub(op(A)) * sub(op(B))) + beta * C
	void pzgemm_    (char* transA, char* transB, int *m, int *n, int *k,
					 double *alpha, cxdb *A, int *ia, int *ja, int * descA, 
					 cxdb *B, int *ib, int *jb, int * descB, double *beta, 
					 double *C, int * ic, int * jc, int *descC);
	void pcgemm_    (char* transA, char* transB, int *m, int *n, int *k,
					 double *alpha, cxfl *A, int *ia, int *ja, int * descA, 
					 cxfl *B, int *ib, int *jb, int * descB, double *beta, 
					 double *C, int * ic, int * jc, int *descC);

}


// C++ convenience
inline int indxl2g (int lidx, int nb, int iproc, int isrcproc, int nprocs) {
	int fortran_idxloc = lidx + 1;
	return indxl2g_ (&fortran_idxloc, &nb, &iproc, &isrcproc, &nprocs) - 1;
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
	gemm  (char* transA, char* transB, int *m, int *n, int *k,	double *alpha, 
		   cxfl *A, int *ia, int *ja, int *descA, cxfl *B, int *ib, int *jb, 
		   int *descB, double *beta, double *C, int * ic, int * jc, int *descC) {
		pcgemm_  (transA, transB, m, n, k, alpha, A, ia, ja, descA, B, ib, jb, 
				  descB, beta, C, ic, jc, descC);
	}

    inline static void 
	gesvd (char* jbu, char* jbvt, int* m, int* n, cxfl* A, int* ia, int* ja, 
		   int* descA, double* s, cxfl* U, int* iu, int* ju, int* descU, 
		   cxfl* VT, int* ivt, int* jvt, int* descVT, cxfl* WORK, int* lwork, 
		   double* rwork, int* info) {
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
	gemm  (char* transA, char* transB, int *m, int *n, int *k,	double *alpha, 
		   cxdb *A, int *ia, int *ja, int *descA, cxdb *B, int *ib, int *jb, 
		   int *descB, double *beta, double *C, int * ic, int * jc, int *descC) {
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
