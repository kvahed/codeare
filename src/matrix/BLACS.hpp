#ifndef __BLACS_HPP__
#define __BLACS_HPP__

#include "config.h"

#ifdef HAVE_MPI

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
    
};

#endif // HAVE_MPI

#endif // __BLACS_HPP__
