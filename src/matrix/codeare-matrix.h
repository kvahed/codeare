/*
 *  codeare Copyright (C) 2007-2010 
 *                        Kaveh Vahedipour
 *                        Forschungszentrum Juelich, Germany
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

#ifndef __MATRIX_H__
#define __MATRIX_H__

#ifdef PARC_MODULE_NAME

  #include     "MrServers/MrVista/include/Ice/IceBasic/IceAs.h"
  #include     "MrServers/MrVista/include/Ice/IceBasic/IceObj.h"
  #include     "MrServers/MrVista/include/Parc/Trace/IceTrace.h"

#else
enum IceDim {
    COL, LIN, CHA, SET, ECO, PHS, REP, SEG, PAR, SLC, IDA, IDB, IDC, IDD, IDE, AVE
};
#endif

#ifndef INVALID_DIM
    #define INVALID_DIM 16
#endif

#include "config.h"
#include "common.h"
#include "OMP.hpp"
#include "Complex.hpp"
#include "SIMD.hpp"

#include "cycle.h"            // FFTW cycle implementation

#include <assert.h>
#include <iostream>
#include <memory>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <valarray>
#include <vector>
#include <ostream>
#include <string>

#include <sys/mman.h>
#include <sys/types.h>


#define ICE_SHRT_MAX 4095


#ifdef HAVE_MAT_H
#include "mat.h"
#endif

#include "Ptr.h"

#define HAVE_VALARRAY 1


/**
 * Short test if the matrix is a vector.
 */
# define VECT(M) assert((M)->width() == 1 || (M)->height() == 1);

 
/**
 * Return minimum of two numbers
 */
# define MAX(A,B) (A > B ? A : B)


/**
 * Return maximum of two numbers 
 */
# define MIN(A,B) (A < B ? A : B)


/**
 * Return rounded value as in MATLAB.
 * i.e. ROUND( 0.1) ->  0
 *      ROUND(-0.5) -> -1
 *      ROUND( 0.5) ->  1
 */
# define ROUND(A) ( floor(A) + ((A - floor(A) >= 0.5) ? (A>0 ? 1 : 0) : 0))


/**
 * @brief   Memory paradigm (share, opencl or message passing)
 */
enum    paradigm {

  SHM, /**< @brief Shared memory (Local RAM) */
  OCL, /**< @brief Open CL GPU RAM */
  MPI  /**< @brief Distributed memory */
  
};

/**
 * @brief   Matrix template.<br/>
 *          Core data structure
 *          
 * @author  Kaveh Vahedipour
 * @date    Mar 2010
 */


template <class T, paradigm P = SHM>
class Matrix  : public SmartObject {
    

public:
    
    
    Matrix ();
    Matrix (const std::vector<size_t>& dim);
	virtual ~Matrix();
    Matrix (const size_t dim[INVALID_DIM]);
    Matrix (const std::vector<size_t>& dim, const std::vector<float>& res);
    Matrix (const size_t dim[INVALID_DIM], const float res[INVALID_DIM]);
    Matrix (const size_t& n);
    Matrix (const size_t& m, const size_t& n);
    Matrix (const size_t& m, const size_t& n, const size_t& k);
    Matrix (const size_t col, 
            const size_t lin, 
            const size_t cha,
            const size_t set,
            const size_t eco = 1,
            const size_t phs = 1,
            const size_t rep = 1,
            const size_t seg = 1,
            const size_t par = 1,
            const size_t slc = 1,
            const size_t ida = 1,
            const size_t idb = 1,
            const size_t idc = 1,
            const size_t idd = 1,
            const size_t ide = 1,
            const size_t ave = 1);
    Matrix (const Matrix<T,P> &M);
    T operator[] (const size_t& p) const;
    T& operator[] (const size_t& p);
    const T* Memory (const size_t p = 0) const;
    inline std::valarray<T>&            
    Container           ();
    inline std::valarray<T>            
    Container           () const;
    inline T            
    At                  (const size_t& p) const;
    inline T&           
    At                  (const size_t& pos) ;
    inline T            
    At                  (const size_t& x, const size_t& y) const;
    virtual const char* GetClassName() const;

};

/*
  

#ifdef PARC_MODULE_NAME

template <class T, paradigm P> long 
Matrix<T,P>::Import     (const IceAs* ias, const size_t pos) {
    
    ICE_SET_FN("Matrix<T,P>::Import(IceAs, long)")
        
    int  i    = 0;
    long size = 1;
    
    for (i = 0; i < INVALID_DIM; i++)
        size *= (ias->getLen(IceDim(i)) <= 1) ? 1 : ias->getLen(IceDim(i));
    
    T* data = (T*) ias->calcSplObjStartAddr() ;
    
    for (i = 0; i < size; i++, data++)
        _M[i+pos] = *data;
    
    return size;
    
}


template <class T, paradigm P> long 
Matrix<T,P>::Import(const IceAs* ias) {
    
    ICE_SET_FN("Matrix<T,P>::Import(IceAs)")
        
    int i;
    
    for (i = 0; i < INVALID_DIM; i++)
        _dim[i] = (ias->getLen(IceDim(i)) <= 1) ? 1 : ias->getLen(IceDim(i));
    
    _M = new T[Size()]();
    nb_alloc = 1;
    
    T* data = (T*) ias->calcSplObjStartAddr() ;
    
    for (i = 0; i < Size(); i++, data++)
        _M[i] = *data;
    
    return Size();
    
}


template <class T, paradigm P> long 
Matrix<T,P>::Export (IceAs* ias) const {
    
    ICE_SET_FN("Matrix<T,P>::Export(IceAs)")
        
    T* data = (T*) ias->calcSplObjStartAddr() ;
    
    for (int i = 0; i < Size(); i++, data++)
        *data = _M[i];
    
    return Size();
    
}


template <class T, paradigm P> long
Matrix<T,P>::Export (IceAs* ias, const size_t pos) const {

    ICE_SET_FN("Matrix<T,P>::Export(IceAs, long)")
        
        int  i    = 0;
    long size = 1;
    
    for (i = 0; i < INVALID_DIM; i++) {
        size *= (ias->getLen(IceDim(i)) <= 1) ? 1 : ias->getLen(IceDim(i));
    }
    
    T* data = (T*) ias->calcSplObjStartAddr() ;
    
    for (i = 0; i < size; i++, data++)
        *data = _M[i+pos];
    
    return size;
    
}

#endif // ICE

#ifdef HAVE_MPI
#include "Grid.hpp"

template<> inline
Matrix<float,MPI>::Matrix (const size_t& cols, const size_t& rows) {

    // Get at home
    int info;
    int izero = 0;
    Grid& gd = *Grid::Instance();
    
    // Global size
    _bs      = 8;
    _gdim[0] = cols;
    _gdim[1] = rows;
		
    // Local size (only with MPI different from global)
    _dim[0] = numroc_ (&_gdim[0], &_bs, &gd.mc, &izero, &gd.nc);
    _dim[1] = numroc_ (&_gdim[1], &_bs, &gd.mr, &izero, &gd.nr);
	
    // Allocate
    this->_M.resize(Size());
	
    // RAM descriptor 
    int dims[2]; dims[0] = _dim[0]; dims[1] = _dim[1];
    descinit_(_desc, &_gdim[0], &_gdim[1], &_bs, &_bs, &izero,
              &izero, &gd.ct, dims, &info);

#ifdef BLACS_DEBUG
    printf ("info(%d) desc({%d, %d, %4d, %4d, %d, %d, %d, %d, %4d})\n", 
            info,     _desc[0], _desc[1], _desc[2], _desc[3], 
            _desc[4], _desc[5], _desc[6], _desc[7], _desc[8]);
#endif

    char a = 'a';
    Cblacs_barrier (gd.ct, &a);

}
#endif
*/

#endif // __MATRIX_H__
