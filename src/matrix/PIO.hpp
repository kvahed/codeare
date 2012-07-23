#ifndef __PIO_HPP__
#define __PIO_HPP__

#include "PMatrix.hpp"

/*
extern "C" {
    void   Cblacs_pinfo    (int* mypnum, int* nprocs);
    void   Cblacs_get      (int context, int request, int* value);
    int    Cblacs_gridinit (int* context, char * order, int np_row, int np_col);
    void   Cblacs_gridinfo (int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
    void   Cblacs_gridexit (int context);
    void   Cblacs_exit     (int error_code);
    void   Cblacs_barrier  (int context, char *scope);
    void   Cigebs2d        (int context, char *scope, char *top, int m, int n, int    *A, int lda);
    void   Cigebr2d        (int context, char *scope, char *top, int m, int n, int    *A, int lda, int rsrc, int csrc);
    void   Cdgebs2d        (int context, char *scope, char *top, int m, int n, double *A, int lda);
    void   Cdgebr2d        (int context, char *scope, char *top, int m, int n, double *A, int lda, int rsrc, int csrc);
    void   Cpdgemr2d       (int M, int N,
                            double *A, int IA, int JA, int *ADESC,
                            double *B, int IB, int JB, int *BDESC,
                            int CTXT);
    int    numroc_         (int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
    void   descinit_       (int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);
    int    descset_        (int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld);
    double pdlamch_        (int *ictxt , char *cmach);
    double pdlange_        (char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);
    void   pdlacpy_        (char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);
    void   pdgesv_         (int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv, double *B, int *ib, int *jb, int *descb, int *info);
    void   pdgemm_         (char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * ALPHA,
                            double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB,
                            double * BETA, double * C, int * IC, int * JC, int * DESCC );
    int    indxg2p_        (int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
    int    indxl2g_        (int *indxloc , int *nb, int *iproc, int *isrcproc, int *nprocs);
    void   pdgeadd_        (char *TRANS, int * M, int * N,
                            double * ALPHA,
                            double * A, int * IA, int * JA, int * DESCA,
                            double * BETA,
                            double * C, int * IC, int * JC, int * DESCC);
} // }}}

inline int IndxL2G(int idxloc, int nb, int iproc, int isrcproc, int nprocs)
{ // {{{

    // (FORTRAN) indxl2g_ function uses indexes which starts at (one)
    int fortran_idxloc = idxloc + 1;
    return indxl2g_(&fortran_idxloc, &nb, &iproc, &isrcproc, &nprocs) - 1;

} // }}}

inline void PDLAFill(LinAlg::Matrix<double> const & AA, double *A, int *DescA, int iRread, int iCread, double *Work)
{ // {{{

   
    // Parameters
    //int block_cyclic_2d = 1;
    //int dtype_          = 1;
    int ctxt_           = 2;
    //int m_              = 3;
    //int n_              = 4;
    int mb_             = 5;
    //int nb_             = 6;
    //int rsrc_           = 7;
    //int csrc_           = 8;
    //int lld_            = 9;
    int one = 1;

    // Local scalars
    bool isioprocessor = false;
    int  ictxt         = DescA[ctxt_-1];
    int  lwork         = DescA[mb_-1];
    int  mycol         = 0;
    int  myrow         = 0;
    int  npcol         = 0;
    int  nprow         = 0;

    // Check if this is the IO processor
    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);
    isioprocessor = ((myrow==iRread)&&(mycol==iCread));

    // Get number of rows and columns
    int iwork[2];
    if (isioprocessor)
    {
        iwork[0] = AA.Rows();
        iwork[1] = AA.Cols();
        //igebs2d_(&ictxt, "All", " ", 2, 1, &iwork, 2);
        Cigebs2d(ictxt, "All", " ", 2, 1, iwork, 2);
    }
    else
    {
        //igebr2d_(&ictxt, "All", " ", 2, 1, &iwork, 2, &iRread, &iCread);
        Cigebr2d(ictxt, "All", " ", 2, 1, iwork, 2, iRread, iCread);
    }
    int m = iwork[0];
    int n = iwork[1];

    // DESCSET initializes a descriptor vector 
    int descwork[9];
    int mm   = Util::Max(1, Util::Min(m, lwork));
    int nn   = Util::Max(1, static_cast<int>(lwork/mm));
    int mb   = mm;
    int nb   = nn;
    int rsrc = iRread;
    int csrc = iCread;
    int ldd  = Util::Max(1, mm);
    descset_(descwork, &mm, &nn, &mb, &nb, &rsrc, &csrc, &ictxt, &ldd);

    // Fill matrix
    for (int jstart=0; jstart<n; jstart+=nn)
    {
        int jend  = Util::Min(n, jstart+nn);
        int jsize = jend - jstart;
        for (int istart=0; istart<m; istart+=mm)
        {
            int    iend  = Util::Min(m, istart+mm);
            int    isize = iend - istart;
            double alpha = 1.0;
            double beta  = 0.0;
            if (isioprocessor)
            {
                for (int j=0; j<jsize; j++)
                for (int i=0; i<isize; i++)
                    Work[i+j*ldd] = AA(i,j);
            }
            pdgeadd_("N", &isize, &jsize,
                     &alpha,
                     Work, &one, &one, descwork,
                     &beta,
                     A, &istart, &jstart, DescA);
        }
    }

    // Flag (?)
    Work[0] = DescA[mb_];

} // }}}
*/

template <class T> inline static void
print (const PMatrix<T>& PM, const std::string& name = "pmat", int nout = 6) {
    
	int   m     = PM.g_m();
	int   n     = PM.g_n();
	T*    iwork = (T*) malloc (PM.m() * sizeof(T));
	int len     = name.length();

	ScalapackTraits<T>::pxlaprnt (&m, &n, PM.Data(), &izero, &izero, PM.Desc(), &izero, &izero, name.c_str(), &nout, iwork, len);

	free (iwork);
    
}

#endif //__PIO_HPP__
