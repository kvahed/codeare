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


#include "linalg/ScalapackTraits.hpp"


#include "config.h"
#include "common.h"
#include "OMP.hpp"
#include "Complex.hpp"
#include "SSE.hpp"

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
    
    
    /**
     * @name Constructors and destructors
     *       Constructors and destructors
     */
    //@{
    
    
    /**
     * @brief           Contruct 1-dim with single element.
     */
    Matrix              ();
    
    
    /**
     * @brief           Construct matrix with dimension array
     *
     * @param  dim      All 16 Dimensions
     */
    Matrix              (const std::vector<size_t>& dim);
    
    /**
     * @brief           Construct matrix with dimension and resolution arrays
     *
     * @param  dim      All 16 Dimensions
     */
    Matrix              (const size_t dim[INVALID_DIM]);
    
    
    /**
     * @brief           Construct matrix with dimension and resolution arrays
     *
     * @param  dim      All 16 Dimensions
     * @param  res      All 16 Resolutions
     */
    Matrix              (const std::vector<size_t>& dim, const std::vector<float>& res);
    
    
    /**
     * @brief           Construct matrix with dimension and resolution arrays
     *
     * @param  dim      All 16 Dimensions
     * @param  res      All 16 Resolutions
     */
    Matrix              (const size_t dim[INVALID_DIM], const float res[INVALID_DIM]);
    
    
    /**
     * @brief           Construct square 2D matrix
     *
     * Usage:
     * @code{.cpp}
     *   Matrix<double> m (5); // Empty 5x5 matrix
     * @endcode
     *
     * @param  n        Rows & Columns
     */
    Matrix              (const size_t& n) ;
    
    
    /**
     * @brief           Construct 2D matrix
     *
     * Usage:
     * @code{.cpp}
     *   Matrix<double> m (5,4); // Empty 5x4 matrix
     * @endcode
     *
     * @param  m        Rows
     * @param  n        Columns
     */
    Matrix              (const size_t& m, const size_t& n);
    
    
    /**
     * @brief           Construct 3D volume
     *
     * Usage:
     * @code{.cpp}
     *   Matrix<double> m (5,4,6); // Empty 5x4x6 matrix
     * @endcode
     *
     * @param  m        Rows
     * @param  n        Columns
     * @param  k        Slices
     */
    Matrix              (const size_t& m, const size_t& n, const size_t& k);
    

    /**
     * @brief           Construct 4D volume
     *
     * @param  col      Scan
     * @param  lin      Phase encoding lines
     * @param  cha      Channels
     * @param  set      Sets
     * @param  eco      Echoes
     * @param  phs      Phases
     * @param  rep      Repetitions
     * @param  seg      Segments
     * @param  par      Partitions
     * @param  slc      Slices
     * @param  ida      IDA
     * @param  idb      IDB
     * @param  idc      IDC
     * @param  idd      IDD
     * @param  ide      IDE
     * @param  ave      Averages
     */
    Matrix              (const size_t col, 
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
    

    
    /**
     * @brief           Delete array containing data.
     */
    virtual
    ~Matrix             ();


    /**
     * @brief           Copy constructor
     *
     * Usage:
     * @code{.cpp}
     *   Matrix<cxfl> m (n); // Copy n into m
     * @endcode
     *
     * @param  M        Right hand side
     */
    Matrix              (const Matrix<T,P> &M);


    //@}



    /**
     * @name            Import export functions for ICE access specifiers.<br/>
     *                  Ice access specifiers can be handled in one of the following ways.<br/>
     *                  It is crucial to understand that the
     */

    //{@

    // Only if compiled within IDEA we know of access specifiers.
#ifdef PARC_MODULE_NAME
    

    /**
     * @brief           Reset and import data from IceAs
     *                   
     * @param  ias      IceAs containing data
     * @return          Amount of data read
     */
    size_t       
    Import              (const IceAs* ias);


    /**
     * @brief           Continue import from IceAs
     *                   
     * @param  ias      IceAs containing data
     * @param  pos      Import data starting at position pos of own repository
     * @return          Amount of data read
     */
    size_t       
    Import              (const IceAs* ias, const size_t pos);


    /**
     * @brief           Import with MDH
     *                   
     * @param  ias      IceAs containing data
     * @param  mdh      Measurement data header      
     * @return          Amount of data read
     */
    size_t       
    Import              (const IceAs* ias, sMDH* mdh);


    /**
     * @brief           Export data to ias
     *                   
     * @param  ias      IceAs for data export
     * @return          Amount of data exported
     */
    size_t         
    Export              (IceAs* ias) const;


    /**
     * @brief           Partially export data to ias 
     * 
     * @param  ias      IceAs for data export
     * @param  pos      Export data starting at position pos of our repository
     */
    size_t
    Export              (IceAs* ias, const size_t pos) const;
 
    #endif

    //@}

    
    /**
     * @name            Elementwise access
     *                  Elementwise access
     */
    
    //@{
    
    
    /**
     * @brief           Get copy of p-th element.
     *
     * Usage:
     * @code{.cpp}
     *   Matrix<float> m = rand<float> (10,1);
     *   float f = m[6]; // right hand side
     * @endcode
     *
     * @param  p        Requested position.
     * @return          Value at _M[p].
     */
    T                   
    operator[]          (const size_t& p) const;
    
    
    /**
     * @brief           Get reference to pth element.
     *
     * Usage:
     * @code{.cpp}
     *   Matrix<float> m = rand<float> (10,1);
     *   m[6] = 1.8; // left hand side
     * @endcode
     *
     * @param  p        Requested position.
     * @return          Reference to _M[p].
     */
    T                   
    &operator[]         (const size_t& p);

    
    /**
     * @brief           Get pointer to data starting from p-th (default:0) element.
     *  
     * @param  p        Position
     *
     * @return          Data 
     */
    inline const T*            
    Memory             (const size_t p = 0)  const {
        assert (p < Size());
        return &(_M[p]);
    }

    
    /**
     * @brief           Get data (lhs)
     *  
     * @return          Data 
     */
    inline std::valarray<T>&            
    Container           ()  {
        return _M;
    }

    
    /**
     * @brief           Get data (rhs)
     *  
     * @return          Data 
     */
    inline std::valarray<T>            
    Container           ()  const {
        return _M;
    }

    
    /**
     * @brief           Element at position p (rhs)
     *  
     * @param  p        Position
     * @return          Value at _M[p]
     */
    inline T            
    At                  (const size_t& p) const {

        assert (p < Size());
        return _M[p];

    }


    /**
     * @brief            Element at position (lhs)
     *  
     * @param  pos       Position
     * @return           Reference to _M[p]
     */
    inline T&           
    At                  (const size_t& pos) {

        assert (pos < Size());
        return _M[pos];

    }


    
    /**
     * @brief           Get element in (first) slice
     *  
     * @param  x        Column
     * @param  y        Line
     *
     * @return          Value
     */
    inline T            
    At                  (const size_t& x, const size_t& y) const {

        assert (x < _dim[0]);
        assert (y < _dim[1]);

        return _M[x + _dim[COL]*y];

    }

    

    /**
     * @brief            Reference to value in slice
     *  
     * @param  x         Column
     * @param  y         Line
     *
     * @return           Reference
     */
    inline T&           
    At                  (const size_t& x, const size_t& y) {

        assert (x < _dim[0]);
        assert (y < _dim[1]);

        return _M[x + _dim[COL]*y];

    }

    
    /**
     * @brief          Get value in volume
     *  
     * @param  x       Column
     * @param  y       Line
     * @param  z       Slice
     *
     * @return         Value
     */
    inline T            
    At                   (const size_t& x, const size_t& y, const size_t& z) const {

    	assert (x < _dim[0]);
    	assert (y < _dim[1]);
    	assert (z < _dim[2]);

        return _M[x + _dim[0]*y + _dim[0]*_dim[1]*z];

    }
    
    

    /**
     * @brief            Reference to value in volume
     *  
     * @param  x         Column
     * @param  y         Line
     * @param  z         Slice
     *
     * @return           Reference
     */
    inline T&            
    At                   (const size_t& x, const size_t& y, const size_t& z) {

    	assert (x < _dim[0]);
    	assert (y < _dim[1]);
    	assert (z < _dim[2]);

        return _M[x + _dim[0]*y + _dim[0]*_dim[1]*z];

    }
    
    

    /**
     * @brief            Get value in store
     *  
     * @param  col       Column
     * @param  lin       Line
     * @param  cha       Channel
     * @param  set       Set
     * @param  eco       Echo
     * @param  phs       Phase
     * @param  rep       Repetition
     * @param  seg       Segment
     * @param  par       Partition
     * @param  slc       Slice
     * @param  ida       Free index A
     * @param  idb       Free index B
     * @param  idc       Free index C
     * @param  idd       Free index D
     * @param  ide       Free index E
     * @param  ave       Average
     * @return           Value at position
     */
    inline T            
    At                   (const size_t& col, 
                          const size_t& lin, 
                          const size_t& cha,
                          const size_t& set,
                          const size_t& eco,
                          const size_t& phs = 0,
                          const size_t& rep = 0,
                          const size_t& seg = 0,
                          const size_t& par = 0,
                          const size_t& slc = 0,
                          const size_t& ida = 0,
                          const size_t& idb = 0,
                          const size_t& idc = 0,
                          const size_t& idd = 0,
                          const size_t& ide = 0,
                          const size_t& ave = 0) const {

    	assert (col < _dim[COL]);
    	assert (lin < _dim[LIN]);
    	assert (cha < _dim[CHA]);
    	assert (set < _dim[SET]);
    	assert (eco < _dim[ECO]);
    	assert (phs < _dim[PHS]);
    	assert (rep < _dim[REP]);
    	assert (seg < _dim[SEG]);
    	assert (par < _dim[PAR]);
    	assert (slc < _dim[SLC]);
    	assert (ida < _dim[IDA]);
    	assert (idb < _dim[IDB]);
    	assert (idc < _dim[IDC]);
    	assert (idd < _dim[IDD]);
    	assert (ide < _dim[IDE]);

        return _M [col+
                   lin*_dim[0]+
                   cha*_dim[0]*_dim[1]+
                   set*_dim[0]*_dim[1]*_dim[2]+
                   eco*_dim[0]*_dim[1]*_dim[2]*_dim[3]+
                   phs*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]+
                   rep*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]+
                   seg*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]+
                   par*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]+
                   slc*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]+
                   ida*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]+
                   idb*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]+
                   idc*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]+
                   idd*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]+
                   ide*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]*_dim[13]+
                   ave*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]*_dim[13]*_dim[14]
                  ];

    }
    
    
    /**
     * @brief            Reference to value in volume
     *  
     * @param  col       Column
     * @param  lin       Line
     * @param  cha       Channel
     * @param  set       Set
     * @param  eco       Echo
     * @param  phs       Phase
     * @param  rep       Repetition
     * @param  seg       Segment
     * @param  par       Partition
     * @param  slc       Slice
     * @param  ida       Free index A
     * @param  idb       Free index B
     * @param  idc       Free index C
     * @param  idd       Free index D
     * @param  ide       Free index E
     * @param  ave       Average
     * @return           Reference to position
     */
    inline T&            
    At                   (const size_t& col, 
                          const size_t& lin, 
                          const size_t& cha,
                          const size_t& set,
                          const size_t& eco = 0,
                          const size_t& phs = 0,
                          const size_t& rep = 0,
                          const size_t& seg = 0,
                          const size_t& par = 0,
                          const size_t& slc = 0,
                          const size_t& ida = 0,
                          const size_t& idb = 0,
                          const size_t& idc = 0,
                          const size_t& idd = 0,
                          const size_t& ide = 0,
                          const size_t& ave = 0) {

    	assert (col < _dim[COL]);
    	assert (lin < _dim[LIN]);
    	assert (cha < _dim[CHA]);
    	assert (set < _dim[SET]);
    	assert (eco < _dim[ECO]);
    	assert (phs < _dim[PHS]);
    	assert (rep < _dim[REP]);
    	assert (seg < _dim[SEG]);
    	assert (par < _dim[PAR]);
    	assert (slc < _dim[SLC]);
    	assert (ida < _dim[IDA]);
    	assert (idb < _dim[IDB]);
    	assert (idc < _dim[IDC]);
    	assert (idd < _dim[IDD]);
    	assert (ide < _dim[IDE]);

        return _M [col+
                   lin*_dim[0]+
                   cha*_dim[0]*_dim[1]+
                   set*_dim[0]*_dim[1]*_dim[2]+
                   eco*_dim[0]*_dim[1]*_dim[2]*_dim[3]+
                   phs*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]+
                   rep*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]+
                   seg*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]+
                   par*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]+
                   slc*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]+
                   ida*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]+
                   idb*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]+
                   idc*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]+
                   idd*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]+
                   ide*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]*_dim[13]+
                   ave*_dim[0]*_dim[1]*_dim[2]*_dim[3]*_dim[4]*_dim[5]*_dim[6]*_dim[7]*_dim[8]*_dim[9]*_dim[10]*_dim[11]*_dim[12]*_dim[13]*_dim[14]
                  ];

    }
    

    /**
     * @brief          Cast operator
     *
     * @return         Cast if possible
     */
    template<class S> operator Matrix<S,P> () const;


    /**
     * @brief           Get the element at position p of the vector, i.e. this(p).
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    T                  
    operator()          (const size_t& p) const {

        return this->At(p);

    }

    
    /**
     * @brief           Get value of pth element of repository.
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    T&                 
    operator()          (const size_t& p) {

        return this->At(p);

    }

    
    /**
     * @brief           Get value in slice
     *
     * @param  x        Column
     * @param  y        Line
     * @return          Value
     */
    inline T 
    operator()          (const size_t& x, const size_t& y) const {

        return this->At(x,y);

    }
    

    /**
     * @brief           Reference to value in slice
     *
     * @param  x        Column
     * @param  y        Line
     * @return          Reference
     */
    inline T&                  
    operator()           (const size_t& x, const size_t& y) {

        return this->At(x,y);

    }
    
    
    /**
     * @brief            Get value in volume
     *
     * @param  x         Column
     * @param  y         Line
     * @param  z         Slice
     *
     * @return           Value
     */
    inline T                  
    operator()           (const size_t& x, const size_t& y, const size_t& z) const {

        return this->At(x,y,z);

    }
    
    
    /**
     * @brief            Reference to value in volume
     *
     * @param  x         Column
     * @param  y         Line
     * @param  z         Slice
     *
     * @return           Reference to _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc]
     */
    inline T&                 
    operator()           (const size_t& x, const size_t& y, const size_t& z) {

        return this->At(x,y,z);

    }
    
    
    /** 
     * @brief            Reference to element
     *
     * @param  col      Column
     * @param  lin      Line
     * @param  cha      Channel
     * @param  set      Set
     * @param  eco      Echo
     * @param  phs      Phase
     * @param  rep      Repetition
     * @param  seg      Segment
     * @param  par      Partition
     * @param  slc      Slice
     * @param  ida      Free index A
     * @param  idb      Free index B
     * @param  idc      Free index C
     * @param  idd      Free index D
     * @param  ide      Free index E
     * @param  ave      Average
     * @return          Reference to position
     */
    inline T&                 
    operator()           (const size_t& col, 
                          const size_t& lin, 
                          const size_t& cha,
                          const size_t& set,
                          const size_t& eco = 0,
                          const size_t& phs = 0,
                          const size_t& rep = 0,
                          const size_t& seg = 0,
                          const size_t& par = 0,
                          const size_t& slc = 0,
                          const size_t& ida = 0,
                          const size_t& idb = 0,
                          const size_t& idc = 0,
                          const size_t& idd = 0,
                          const size_t& ide = 0,
                          const size_t& ave = 0) { 

        return this->At (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave);
        
    }

    /** 
     * @brief            Value at ...
     *
     * @param  col      Column
     * @param  lin      Line
     * @param  cha      Channel
     * @param  set      Set
     * @param  eco      Echo
     * @param  phs      Phase
     * @param  rep      Repetition
     * @param  seg      Segment
     * @param  par      Partition
     * @param  slc      Slice
     * @param  ida      Free index A
     * @param  idb      Free index B
     * @param  idc      Free index C
     * @param  idd      Free index D
     * @param  ide      Free index E
     * @param  ave      Average
     * @return          Value
     */
    inline T
    operator()           (const size_t& col, 
                          const size_t& lin, 
                          const size_t& cha,
                          const size_t& set,
                          const size_t& eco = 0,
                          const size_t& phs = 0,
                          const size_t& rep = 0,
                          const size_t& seg = 0,
                          const size_t& par = 0,
                          const size_t& slc = 0,
                          const size_t& ida = 0,
                          const size_t& idb = 0,
                          const size_t& idc = 0,
                          const size_t& idd = 0,
                          const size_t& ide = 0,
                          const size_t& ave = 0) const { 
        
        return this->At (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave);

    }

    //@}
    


    /**
     * @name            Friend operators
     *                  Who doesn't need friends
     */
    
    //@{
    

    //--
    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator*  (const double& s, const Matrix<T,P>& m) { 
        return   m * s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator*  (const float& s, const Matrix<T,P> &m) { 
        return   m * s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator*  (const short& s, const Matrix<T,P> &m) { 
        return   m * s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator*  (const long& s, const Matrix<T,P> &m) { 
        return   m * s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator*  (const cxfl& s, const Matrix<T,P> &m) { 
        return   m * s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator*  (const cxdb& s, const Matrix<T,P> &m) { 
        return   m * s;
    }


    //--
    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator+  (const double& s, const Matrix<T,P> &m) { 
        return   m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator+  (const float& s, const Matrix<T,P> &m) { 
        return   m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator+  (const short& s, const Matrix<T,P> &m) { 
        return   m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator+  (const long& s, const Matrix<T,P> &m) { 
        return   m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator+  (const cxfl& s, const Matrix<T,P> &m) { 
        return   m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator+  (const cxdb& s, const Matrix<T,P> &m) { 
        return   m + s;
    }


    //--
    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator-  (const double& s, const Matrix<T,P> &m) { 
        return -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator-  (const float& s, const Matrix<T,P> &m) { 
        return -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator-  (const short& s, const Matrix<T,P> &m) { 
        return -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator-  (const long& s, const Matrix<T,P> &m) { 
        return   -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator-  (const cxfl& s, const Matrix<T,P> &m) { 
        return   -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator-  (const cxdb& s, const Matrix<T,P> &m) { 
        return   -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator/  (const double& s, const Matrix<T,P> &m) { 

    	assert (s != 0.0);

        if (s == 1.0)
            return Matrix<T,P> (m);
        
        Matrix<T,P> res = m;
        res.Container() = s / res.Container();
        return res;
        
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>    
    operator/  (const float& s, const Matrix<T,P> &m) { 

    	assert (s != 0.0);

		if (s == 1.0)
			return Matrix<T,P> (m);

		Matrix<T,P> res = m;
		res.Container() = s / res.Container();
		return res;

    }


    /**
     * @brief           Elementwise equality with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m == s
     */
    inline friend Matrix<bool> 
    operator== (const T& s, const Matrix<T,P>& m) {
        return   m == s;
    }


    /**
     * @brief           Elementwise >= with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m <= t
     */
    inline friend Matrix<bool> 
    operator>= (const T& s, const Matrix<T,P>& m) {
        return   m <= s;
    }


    /**
     * @brief           Elementwise <= with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T<=M
     */
    inline friend Matrix<bool> 
    operator<= (const T& s, const Matrix<T,P>& m) {
        return   m >= s;
    }


    /**
     * @brief           Elementwise unequality with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T!=M
     */
    inline friend Matrix<bool> 
    operator!= (const T& s, const Matrix<T,P>& m) {
        return   m != s;
    }


    /**
     * @brief           Elementwise equality with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T+M
     */
    inline friend Matrix<bool> 
    operator>  (const T& s, const Matrix<T,P>& m) {
        return   m <  s;
    }


    /**
     * @brief           Elementwise < with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T+M
     */
    inline friend Matrix<bool> 
    operator<  (const T& s, const Matrix<T,P>& m) {
        return   m >  s;
    }


    /**
     * @brief           Elementwise equality with scalar (lhs)
     *
     * @param  mb       Scalar lhs
     * @param  m        Matrix rhs
     * @return          T+M
     */
    inline friend Matrix<T,P>    
    operator&  (const Matrix<bool>& mb, const Matrix<T,P>& m) {
        return   m & mb;
    }

    //@}



    //@}
    
    /**
     * @name            Dimensions
     *                  Some convenience functions to access dimensionality
     */
    
    //@{

    /**
     * @brief           Get number of rows, i.e. tmp = size(this); tmp(1).
     *
     * @return          Number of rows.
     */
    inline size_t                 
    Height              () const {
        return _dim[0];
    }
    
    
    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    inline size_t                 
    Width               () const {
        return _dim[1];
    }
    
#ifdef HAVE_MPI    
    /**
     * @brief           Get number of rows, i.e. tmp = size(this); tmp(1).
     *
     * @return          Number of rows.
     */
    inline size_t                 
    GHeight              () const {
        return _gdim[0];
    }
    
    
    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    inline size_t                 
    GWidth               () const {
        return _gdim[1];
    }


    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    inline const int*
    Desc               () const {
        return _desc;
    }
#endif
    

    /**
     * @brief           Get resolution a given dimension.
     *
     * @param   i       Dimension
     * @return          Resolution .
     */
    inline float          
    Res                 (const size_t& i) const {
        assert (i < INVALID_DIM);
        return _res[i];
    }
    
    
    /**
     * @brief           Rresolution a given dimension.
     *
     * @param   i       Dimension
     * @return          Resolution
     */
    inline float&          
    Res                 (const size_t& i)       {
        assert (i < INVALID_DIM);
        return _res[i];
    }
    


    /**
     * @brief           Resolution array
     *
     * @return          All resolutions
     */
    const float*
    Res                 () const {
        return &_res[0];
    }
    

    
    /**
     * @brief           Get size a given dimension.
     *
     * @param   i       Dimension
     * @return          Dimension
     */
    inline size_t          
    Dim                 (const size_t& i)  const {
        assert (i < INVALID_DIM);
        return _dim[i];
    }
    
    
    /**
     * @brief           Get dimension array
     *
     * @return          All dimensions
     */
    inline const size_t*   
    Dim                 ()                  const {
        return _dim;
    }
    

    /**
     * @brief           Get dimension vector
     *
     * @return          All dimensions
     */
    inline std::vector<size_t>   
    DimVector           ()                  const {
		std::vector<size_t> dim (INVALID_DIM,1);
		for (size_t i = 0; i < INVALID_DIM; i++)
			dim[i] = _dim[i];
        return dim;
    }
    

    /**
     * @brief           Get size a given dimension.
     *
     * @return          Number of rows.
     */
    inline size_t          
    Dim                 (const int& i)      const {
        assert (i < INVALID_DIM);
        return _dim[i];
    }
    
    
    /**
     * @brief        Resize to mxn 2D matrix.<br/>Preserving data while shrinking.<br/>Adding zeros when growing.
     *               
     * @param  m     # Rows
     * @param  n     # Columns
     */
    inline void         
    Resize          (const size_t& m, const size_t& n)                               {

        _dim[0] = m;
        _dim[1] = n;

        for (size_t i = 2; i < INVALID_DIM; i++)
            _dim[i] = 1;

		std::valarray<T> tmp = _M;

        _M.resize(Size(), T(0));

        memcpy (&_M[0], &tmp[0], MIN(Size(),tmp.size()) * sizeof(T));

    }
    

    /**
     * @brief           Purge data and free RAM.
     */
    inline void         
    Clear               ()                                      {
        
        for (size_t i = 0; i < INVALID_DIM; i++)
            _dim[i] = 1;

        _M.resize(1);

    }
    
    
    /**
     * @brief           Reset. i.e. Set all fields = T(0)
     */
    inline void         
    Zero               ()                                      {
        
        _M = T(0); 
        
    }
    

    //@}
    
    

    /**
     * @name            Some operators
     *                  Operator definitions. Needs big expansion still.
     */
    
    //@{
    

    
    /**
     * @brief           Assignment operator. i.e. this = m.
     *
     * @param  M        The assigned matrix.
     */
    Matrix<T,P>&
    operator=           (const Matrix<T,P>& M);
    
    
    /**
     * @brief           Assignment operator. i.e. this = m.
     *
     * @param  v        Data vector (size must match numel(M)).
     */
    Matrix<T,P>&
    operator=           (const std::valarray<T>& v);
    
    
    /**
     * @brief           Assignment operator. Sets all elements s.
     *
     * @param  s        The assigned scalar.
     */
    Matrix<T,P>&
    operator=           (const T& s) {
        this->_M = s;
        return *this;
    }
    
    
    /**
     * @brief           Matrix product. i.e. this * M.
     *
     * @param  M        The factor.
     */
    Matrix<T,P>           
    operator->*         (const Matrix<T,P>& M) const;
   
    
    /**
     * @brief           Elementwise substruction of two matrices
     *
     * @param  M        Matrix substruent.
     */
    template <class S> Matrix<T,P>           
    operator-           (const Matrix<S,P>& M) const;
    
    
    /**
     * @brief           Elementwise substruction all elements by a scalar
     *
     * @param  s        Scalar substruent.
     */
    template <class S> Matrix<T,P>           
    operator-           (const S& s) const;
    
    
    /**
     * @brief           ELementwise substraction and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
     * @return          Result
     */
    template <class S>  Matrix<T,P>&           
    operator-=          (const Matrix<S,P>& M);
    
    
    /**
     * @brief           ELementwise substration with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
     * @return          Result
     */
    template <class S> Matrix<T,P>&
    operator-=         (const S& s);
    
    
    /**
     * @brief           Unary minus (additive inverse)
     *
     * @return          Negation
     */
    Matrix<T,P>           
    operator-           () const;
    
    
    /**
     * @brief           Unary plus
     *
     * @return          Identity
     */
    Matrix<T,P>           
    operator+           () const;
    
    
    /**
     * @brief           Elementwise addition of two matrices
     *
     * @param  M        Matrix additive.
     */
    template <class S> Matrix<T,P>           
    operator+          (const Matrix<S,P>& M) const;
    
    
    /**
     * @brief           Elementwise addition iof all elements with a scalar
     *
     * @param  s        Scalar additive.
     */
    template <class S> Matrix<T,P>           
    operator+           (const S& s) const;
    
    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
     * @return          Result
     */
    template <class S> Matrix<T,P>&
    operator+=          (const Matrix<S,P>& M);
    
    
    /**
     * @brief           ELementwise addition with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
     * @return          Result
     */
    template <class S > Matrix<T,P>&           
    operator+=          (const S& s);
    
    
    /**
     * @brief           Transposition / Complex conjugation. i.e. this'.
     *
     * @return          Matrix::tr()
     */
    Matrix<T,P>           
    operator!           () const;
    
    
    /**
     * @brief           Return a matrix with result[i] = (m[i] ? this[i] : 0).
     *
     * @param  M        The operand
     * @return          Cross-section or zero
     */
    Matrix<T,P>           
    operator&           (const Matrix<bool>& M) const ;
    
    
    /**
     * @brief           Scalar equality. result[i] = (this[i] == m).
     *
     * @param  s        Comparing scalar.
     * @return          Matrix of true where elements are equal s and false else.
     */
    Matrix<bool>        
    operator==          (const T& s) const ;
    
    
    /**
     * @brief           Scalar inequality. result[i] = (this[i] != m). i.e. this ~= m
     *
     * @param  s        Comparing scalar.
     * @return          Matrix of false where elements are equal s and true else.
     */
    Matrix<bool>        
    operator!=          (const T& s) const ;
    
    
    /**
     * @brief           Scalar greater comaprison, result[i] = (this[i] > m). i.e. this > m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    Matrix<bool>        
    operator>           (const T& s) const ;
    
    
    /**
     * @brief           Scalar greater or equal comparison. result[i] = (this[i] >= m). i.e. this >= m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    Matrix<bool>        
    operator>=          (const T& s) const;
    
    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] <= m). i.e. this <= m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    Matrix<bool>        
    operator<=          (const T& s) const;
    
    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] < m). i.e. this < m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    Matrix<bool>        
    operator<           (const T& s) const;
    
    
    /**
     * @brief           Elementwise equality, result[i] = (this[i] == m[i]). i.e. this == m
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    Matrix<bool>        
    operator==          (const Matrix<T,P>& M) const;
    
    
    /**
     * @brief           Elementwise equality, result[i] = (this[i] != m[i]). i.e. this ~= m
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    Matrix<bool>        
    operator!=          (const Matrix<T,P>& M) const;
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] >= m[i]). i.e. this >= m
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    Matrix<bool>        
    operator>=          (const Matrix<T,P>& M) const;
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] <= m[i]). i.e. this <= m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    Matrix<bool>        
    operator<=          (const Matrix<T,P>& M) const;
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] > m[i]). i.e. this > m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    Matrix<bool>        
    operator>           (const Matrix<T,P>& M) const;
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] < m[i]). i.e. this < m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    Matrix<bool>        
    operator<           (const Matrix<T,P>& M) const; 
    
    
    /**
     * @brief           Matrix comparison, result[i] = (m[i] || this[i] ? 1 : 0). i.e. this | m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    Matrix<bool>           
    operator||          (const Matrix<T,P>& M) const;
    
    
    /**
     * @brief           Matrix comparison, result[i] = (m[i] && this[i] ? 1 : 0). i.e. this & m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    Matrix<bool>           
    operator&&          (const Matrix<T,P>& M) const;


    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Power.
     * @return          Result
     */
    Matrix<T,P>           
    operator^           (const float& p) const;
    
    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Power.
     * @return          Result
     */
    Matrix<T,P>&
    operator^=          (const float& p);
    

    /**
     * @brief           Elementwise multiplication. i.e. this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    template <class S> Matrix<T,P>           
    operator*          (const Matrix<S,P> &M) const ;


    /**
     * @brief           Elementwise multiplication with a scalar. i.e. this * m.
     *
     * @param  s        Factor scalar
     * @return          Result
     */
    template <class S> Matrix<T,P>           
    operator*          (const S& s) const ;

    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    template <class S> Matrix<T,P>&
    operator*=         (const Matrix<S,P>& M);
    
    
    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this *= s.
     *
     * @param  s        Factor scalar.
     * @return          Result
     */
    template <class S> Matrix<T,P>&
    operator*=         (const S& s);
    
    
    /**
     * @brief           Elelemtwise division by M.
     *
     * @param  M        The divisor.
     * @return          Result
     */
    template <class S>  Matrix<T,P>           
    operator/          (const Matrix<S,P>& M) const;

    
    /**
     * @brief           Elementwise division by scalar. i.e. this * m.
     *
     * @param  s        The divisor.
     * @return          Result
     */
    template <class S> Matrix<T,P>           
    operator/           (const S& s) const;
    
    /**
     * @brief           ELementwise division and assignment operator. i.e. this = this ./ M.
     *
     * @param  M        Divisor matrix.
     * @return          Result
     */
    template <class S> Matrix<T,P>&
    operator/=         (const Matrix<S,P> &M);
    
    
    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Divisor scalar.
     * @return          Result
     */
    template <class S> Matrix<T,P>&
    operator/=         (const S& s);
    
    
    //@}



    /**
     * @name            Other functions.
     *                  Other functions.
     */
    
    //@{
    

    /**
     * @brief           Who are we?
     *
     * @return          Class name
     */ 
    const char* 
    GetClassName        () const { 
        return _name.c_str(); 
    }


    /**
     * @brief           Who are we?
     *
     * @return          Class name
     */ 
    void 
    SetClassName        (const char* name) { 
        _name = name; 
    }


    /**
     * @brief           Matrix Product.
     *
     * @param   M       The factor
     * @param   transa  Transpose ('T') / Conjugate transpose ('C') the left matrix. Default: No transposition'N'
     * @param   transb  Transpose ('T') / Conjugate transpose ('C') the right matrix. Default: No transposition 'N'
     * @return          Product of this and M.
     */
    Matrix<T,P>           
    prod                (const Matrix<T,P> &M, const char& transa = 'N', const char& transb = 'N') const;
    

    /**
     * @brief           Complex conjugate left and multiply with right.
     *
     * @param   M       Factor
     * @return          Product of conj(this) and M.
     */
    Matrix<T,P>           
    prodt               (const Matrix<T,P> &M) const;
    

    /**
     * @brief           Complex conjugate right and multiply with right.
     *
     * @param   M       Factor
     * @return          Product of conj(this) and M.
     */
    Matrix<T,P>           
    tprod               (const Matrix<T,P> &M) const;
    

    /**
     * @brief           Scalar product (complex: conjugate first vector) using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTC and XDOT
     *
     * @param  M        Factor
     * @return          Scalar product
     */
    T
    dotc (const Matrix<T,P>& M) const;
    
    
    /**
     * @brief           Scalar product using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTU and XDOT
     *
     * @param  M        Factor
     * @return          Scalar product
     */
    T
    dotu (const Matrix<T,P>& M) const;
    
    /**
     * @brief           Scalar product using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTU and XDOT
     *
     * @param  M        Factor
     * @return          Scalar product
     */
    T
    dot (const Matrix<T,P>& M) const;
    
    /**
     * @brief           Transposition / Complex conjugation and transposition.
     *
     * @return          The transposed matrix
     */
    Matrix<T,P>           
    tr   ()             const;
    

    /**
     * @brief           Number of elements
     *
     * @return          Size
     */
    virtual inline size_t
    Size() const {
        
        
        long size = 1;
        
        for (size_t i = 0; i < INVALID_DIM; i++)
            size *= _dim[i];
        
        return size;
        
    }

    //@}


protected:
	
    // Structure
    size_t              _dim[INVALID_DIM]; /// Dimensions
    float               _res[INVALID_DIM]; /// Resolutions

	// Data
    std::valarray<T>    _M;

    // Name
    std::string         _name; 
    
#ifdef HAVE_MPI
    // BLACS 
	int              _bs;
	int              _desc[9]; /**< @brief matrix grid vector */
	int              _gdim[2]; /**< @brief Global dimensions */
#endif
    
    /**
     * @brief           Adjust and resize for Syngo read
     *
     * @param  fname    Syngo MR meas file name
     * @return          Success
     */
    bool
    RSAdjust            (const std::string& fname);

    /* Who do we support? */
    void Validate (short&  t) const {};
    void Validate (long&   t) const {};
    void Validate (size_t& t) const {};
    void Validate (float&  t) const {};
    void Validate (double& t) const {};
    void Validate (cxfl&   t) const {};
    void Validate (cxdb&   t) const {};
	void Validate (bool&   t) const {};

};



#include "linalg/Lapack.hpp"

template <class T, paradigm P> Matrix<T,P> 
Matrix<T,P>::prodt (const Matrix<T,P> &M) const {
    
    return gemm (*this, M, 'C');
    
}


template <class T, paradigm P> Matrix<T,P> 
Matrix<T,P>::prod (const Matrix<T,P> &M, const char& transa, const char& transb) const {
    
    return gemm (*this, M, transa, transb);
    
}


template<class T, paradigm P> T 
Matrix<T,P>::dotc (const Matrix<T,P>& M) const {

    return DOTC (*this, M);
    
}


template<class T, paradigm P>  T 
Matrix<T,P>::dot (const Matrix<T,P>& M) const {
    
    return DOT  (*this, M);
    
}


template <class T, paradigm P> inline 
Matrix<T,P>::Matrix () {

    T t;
    Validate (t);

    for (size_t i = 0; i < INVALID_DIM; i++) {
        _dim [i] = 1;
        _res [i] = 1.0;
    }

    _M.resize(Size(), T(0));
        
    _name = "matrix";

}



template <class T, paradigm P> inline 
Matrix<T,P>::Matrix (const size_t& n) {

	assert (n);

    T t;
    Validate (t);

    _dim [COL] = n;
    _dim [LIN] = n;

    for (size_t i = 2; i < INVALID_DIM; i++)
        _dim [i] = 1;

    for (size_t i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;
    
    _M.resize(n*n, T(0));
    
    _name = "matrix";

}



template <class T, paradigm P> inline 
Matrix<T,P>::Matrix (const size_t& m, const size_t& n) {

	assert (m * n);

    T t;
    Validate (t);

    _dim [0] = m;
    _dim [1] = n;

    for (size_t i = 2; i < INVALID_DIM; i++) 
        _dim [i] = 1;
    
    for (size_t i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;

    _M.resize(m*n, T(0));

	_name = "matrix";

}



template <class T, paradigm P> inline 
Matrix<T,P>::Matrix (const size_t& m, const size_t& n, const size_t& k) {

	assert (m * n * k);

    T t;
    Validate (t);

    _dim [0] = m;
    _dim [1] = n;
    _dim [2] = k;
    
    for (size_t i = 3; i < INVALID_DIM; i++)
        _dim [i] = 1;
    
    for (size_t i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;
    
    _M.resize(m*n*k, T(0));
    


}



template <class T, paradigm P> inline 
Matrix<T,P>::Matrix (const size_t col, const size_t lin, const size_t cha, const size_t set, 
                   const size_t eco, const size_t phs, const size_t rep, const size_t seg, 
                   const size_t par, const size_t slc, const size_t ida, const size_t idb, 
                   const size_t idc, const size_t idd, const size_t ide, const size_t ave) {
    
	assert (col * lin * cha * set * eco * phs * rep * seg * 
			par * slc * ida * idb * idc * idd * ide * ave);

    T t;
    Validate (t);

    _dim[COL] = col;
    _dim[LIN] = lin;
    _dim[CHA] = cha;
    _dim[SET] = set;
    _dim[ECO] = eco;
    _dim[PHS] = phs;
    _dim[REP] = rep;
    _dim[SEG] = seg;
    _dim[PAR] = par;
    _dim[SLC] = slc;
    _dim[IDA] = ida;
    _dim[IDB] = idb;
    _dim[IDC] = idc;
    _dim[IDD] = idd;
    _dim[IDE] = ide;
    _dim[AVE] = ave;
    
    for (size_t i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;
    
    _M.resize(Size(), T(0));
    
}


template <class T, paradigm P> inline 
Matrix<T,P>::Matrix (const std::vector<size_t>& dim) {
    
	assert (dim.size() <= INVALID_DIM);

	size_t n = 1, i = 0;
	
	for (; i < dim.size(); i++)
		n *= dim[i];

	assert (n);

    T t;
    Validate (t);

    for (i = 0; i < dim.size(); i++) {
        _dim[i] = dim[i];
        _res[i] = 1.0;
    }

	for (; i < INVALID_DIM; i++) {
        _dim[i] = 1;
        _res[i] = 1.0;
	}
    
    _M.resize(Size(),T(0));
    
}


template <class T, paradigm P> inline 
Matrix<T,P>::Matrix (const std::vector<size_t>& dim, const std::vector<float>& res) {
    
	assert (dim.size() <= INVALID_DIM);
	assert (dim.size() == res.size());

	size_t n = 1, i = 0;
	
	for (; i < dim.size(); i++)
		n *= dim[i];

	assert (n);

    T t;
    Validate (t);

	for (i = 0; i < dim.size(); i++) {
		_dim[i] = dim[i];
		_res[i] = res[i];
	}

	for (; i < INVALID_DIM; i++) {
		_dim[i] = 1;
		_res[i] = 1.0;
	}
    
    _M.resize(Size(),T(0));

}


template <class T, paradigm P> inline 
Matrix<T,P>::Matrix (const size_t dim[INVALID_DIM]) {
    
	size_t n = 1, i = 0;
	
	for (; i < INVALID_DIM; i++)
		n *= dim[i];

	assert (n);

    T t;
    Validate (t);

	for (i = 0; i < INVALID_DIM; i++) {
		_dim[i] = dim[i];
		_res[i] = 1.0;
	}

    _M.resize(Size(),T(0));

}


template <class T, paradigm P> inline 
Matrix<T,P>::Matrix (const size_t dim[INVALID_DIM], const float res[INVALID_DIM]) {
    
	size_t n = 1, i = 0;
	
	for (; i < INVALID_DIM; i++)
		n *= dim[i];

	assert (n);

    T t;
    Validate (t);

	for (i = 0; i < INVALID_DIM; i++) {
		_dim[i] = dim[i];
		_res[i] = res[i];
	}

    _M.resize(Size(),T(0));

}


template <class T, paradigm P> inline 
Matrix<T,P>::Matrix (const Matrix<T,P> &M) {
    
	if (this != &M) { 

		T t;
		Validate (t);
		
		memcpy (_dim, M.Dim(), INVALID_DIM * sizeof(size_t));
		memcpy (_res, M.Res(), INVALID_DIM * sizeof( float));
		
		_M = M.Container();
		
	}

}



template <class T, paradigm P> inline 
Matrix<T,P>::~Matrix() {};



template <class T, paradigm P> inline Matrix<T,P>&
Matrix<T,P>::operator= (const Matrix<T,P>& M) {
    
    if (this != &M) {

        memcpy (_dim, M.Dim(), INVALID_DIM * sizeof(size_t));
        memcpy (_res, M.Res(), INVALID_DIM * sizeof( float));
        
        _M = M.Container();
        
    }

    return *this;
    
}


template <class T, paradigm P> inline Matrix<T,P>&
Matrix<T,P>::operator= (const std::valarray<T>& v) {
    
	assert (_M.size() == v.size());

    if (&_M != &v)
        _M = v;

    return *this;
    
}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator== (const T& s) const {

    Matrix<bool> res(_dim);
    res.Container() = (_M == s);
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator>= (const T& s) const {

    Matrix<bool> res(_dim);
    res.Container() = (_M >= s);
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator<= (const T& s) const {

    Matrix<bool> res(_dim);
    res.Container() = (_M <= s);
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator!= (const T& s) const {

    Matrix<bool> res(_dim);
    res.Container() = (_M != s);
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator< (const T& s) const {

    Matrix<bool> res(_dim);
    res.Container() = (_M < s);
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator> (const T& s) const {

    Matrix<bool> res(_dim);
    res.Container() = (_M > s);
    return res;

}


template <class T, paradigm P> inline Matrix<T,P> 
Matrix<T,P>::operator->* (const Matrix<T,P> &M) const {

    return this->prod(M);

}


template <class T, paradigm P> inline Matrix<T,P> 
Matrix<T,P>::operator!() const {
	
    for (size_t i = 2; i < INVALID_DIM; i++)
        assert (_dim[i] == 1);
	
    Matrix<T,P> res (_dim[1],_dim[0]);

#pragma omp parallel
	{
#pragma omp for
		for (size_t i = 0; i < _dim[1]; i++)
			for (size_t j = 0; j < _dim[0]; j++)
				res(i,j) = this->At(j,i);
	}
	
    return res;
	
}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator&& (const Matrix<T,P>& M) const {

    Matrix<bool> res(_dim);
    res.Container() = (_M && M.Container());
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator|| (const Matrix<T,P>& M) const {

    for (size_t i=0; i < INVALID_DIM; i++)
        assert (Dim(i) == M.Dim(i));

    Matrix<bool> res(_dim);
    res.Container() = (_M || M.Container());
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator== (const Matrix<T,P>& M) const {

    for (size_t i=0; i < INVALID_DIM; i++)
        assert (Dim(i) == M.Dim(i));

    Matrix<bool> res(_dim);
    res.Container() = (_M == M.Container());
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator>= (const Matrix<T,P>& M) const {

    for (size_t i=0; i < INVALID_DIM; i++)
        assert (Dim(i) == M.Dim(i));

    Matrix<bool> res(_dim);
    res.Container() = (_M >= M.Container());
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator<= (const Matrix<T,P>& M) const {

    for (size_t i=0; i < INVALID_DIM; i++)
        assert (Dim(i) == M.Dim(i));

    Matrix<bool> res(_dim);
    res.Container() = (_M <= M.Container());
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator!= (const Matrix<T,P>& M) const {

    for (size_t i=0; i < INVALID_DIM; i++)
        assert (Dim(i) == M.Dim(i));

    Matrix<bool> res(_dim,_res);
    res.Container() = (_M != M.Container());
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator> (const Matrix<T,P>& M) const {
    
    for (size_t i=0; i < INVALID_DIM; i++)
        assert (Dim(i) == M.Dim(i));

    Matrix<bool> res(_dim,_res);
    res.Container() = (_M > M.Container());
    return res;

}


template <class T, paradigm P> inline Matrix<bool> 
Matrix<T,P>::operator< (const Matrix<T,P>& M) const {

    for (size_t i=0; i < INVALID_DIM; i++)
        assert (Dim(i) == M.Dim(i));

    Matrix<bool> res(_dim,_res);
    res.Container() = (_M < M.Container());
    return res;

}


template <class T, paradigm P> inline Matrix<T,P> 
Matrix<T,P>::operator- () const {

    Matrix<T,P> res (_dim,_res);
    res.Container() = -_M;
    return res;

}


template <class T, paradigm P> inline Matrix<T,P>
Matrix<T,P>::operator+ () const {

    return *this;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P> 
Matrix<T,P>::operator- (const Matrix<S,P> &M) const {

    for (size_t i=0; i < INVALID_DIM; i++)
        assert (Dim(i) == M.Dim(i));

    Matrix<T,P> res = *this;
	res.Container() -= M.Container();
    return res;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P> 
Matrix<T,P>::operator- (const S& s) const {

    Matrix<T,P> res = *this;
	res.Container() -= T(s);
    return res;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P> 
Matrix<T,P>::operator+ (const Matrix<S,P> &M) const {

    for (size_t i=0; i < INVALID_DIM; i++)
        assert (Dim(i) == M.Dim(i));

    Matrix<T,P> res = M;

#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			res[i] += _M[i];
	}

    return res;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P> 
Matrix<T,P>::operator+ (const S& s) const {

    Matrix<T,P> res = *this;
	T t = T(s);

#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			res[i] += t;
	}

    return res;

}


template <class T, paradigm P> inline Matrix<T,P> 
Matrix<T,P>::operator^ (const float& p) const {
    
	Matrix<T,P> res = *this;

#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			res[i] = (p == 0) ? T(1) : pow(res[i],  p);
	}

    return res;

}


template <class T, paradigm P> inline Matrix<T,P>&
Matrix<T,P>::operator ^= (const float& p) {

    size_t i = Size();

#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			_M[i] = (p == 0) ? T(1) : pow(_M[i],  p);
	}

    return *this;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P> 
Matrix<T,P>::operator* (const Matrix<S,P> &M) const {

    for (size_t i = 0; i < INVALID_DIM; i++)
        assert (_dim[i] == M.Dim(i));

    Matrix<T,P> res = M;

#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			res[i] *= _M[i];
	}

	//SSE::process<T,P> (&_M[0], &res[0], M.Size(), SSE::mul<T,P>(), &res[0]) ;

    return res;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P> 
Matrix<T,P>::operator* (const S& s) const {
    
    Matrix<T,P> res = *this; 
	T t = T(s);

#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			res[i] *= t;
	}

    return res;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P>&
Matrix<T,P>::operator *= (const Matrix<S,P> &M) {
    
    size_t i;

    for (i = 0; i < INVALID_DIM; i++)
        assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			_M[i] *= M[i];
	}

	//SSE::process<T,P> (&_M[0], &tmp[0], M.Size(), SSE::mul<T,P>(), &_M[0]);
	
    return *this;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P>&
Matrix<T,P>::operator *= (const S& s) {
    
	T t = T (s);
    
#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			_M[i] *= t;
	}

    return *this;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P>&
Matrix<T,P>::operator += (const Matrix<S,P> &M) {

    for (size_t i = 0; i < INVALID_DIM; i++)
        assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			_M[i] += M[i];
	}

	//SSE::process<T,P> (&_M[0], &tmp[0], M.Size(), SSE::add<T,P>(), &_M[0]);

    return *this;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P>&
Matrix<T,P>::operator+= (const S& s) {

	T t = T (s);
    
#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			_M[i] += t;
	}

    return *this;

}

template <class T, paradigm P> template <class S> inline Matrix<T,P>&
Matrix<T,P>::operator-= (const Matrix<S,P>& M) {

    for (size_t i = 0; i < INVALID_DIM; i++)
        assert (_dim[i] == M.Dim(i));
	
#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			_M[i] -= M[i];
	}

    return *this;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P>&
Matrix<T,P>::operator-= (const S& s) {
    
#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			_M[i] -= T(s);
	}

    return *this;

}


template<class T, paradigm P> template<class S> inline Matrix<T,P>&
Matrix<T,P>::operator /= (const Matrix<S,P> &M) {
    
    size_t i;

    for (i = 0; i < INVALID_DIM; i++)
        assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			_M[i] /= M[i];
	}

    return *this;

}


template <class T, paradigm P> template<class S> inline Matrix<T,P>&
Matrix<T,P>::operator/= (const S& s) {
    
	T zero = T(0.0);
    assert (s != zero);

#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			_M[i] /= T(s);
	}

    return *this;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P> 
Matrix<T,P>::operator/ (const Matrix<S,P>& M) const {

    size_t i;

    for (i = 0; i < INVALID_DIM; i++)
        assert (_dim[i] == M.Dim(i));

    Matrix<T,P> res = *this;
	
#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			res[i] = (M[i] != (T)0) ? _M[i] / M[i] : 0;
	}

    return res;

}


template <class T, paradigm P> template <class S> inline Matrix<T,P> 
Matrix<T,P>::operator/ (const S& s) const {
    
    assert (cabs(s) != 0.0);
	T t = T (s);
    Matrix<T,P> res = *this;
#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			res[i] /= t;
	}

    return res;

}

template <class T, paradigm P> inline T 
Matrix<T,P>::operator[]  (const size_t& p) const {
    
    assert(p <  Size());
    return _M[p];
    
}


template <class T, paradigm P> inline T&
Matrix<T,P>::operator[] (const size_t& p) {
    
    assert(p <  Size());
    return _M[p];
    
}


template<class T, paradigm P> template<class S> inline
Matrix<T,P>::operator Matrix<S,P> () const {

    Matrix<S,P> m (_dim[ 0],_dim[ 1],_dim[ 2],_dim[ 3],
				 _dim[ 4],_dim[ 5],_dim[ 6],_dim[ 7],
				 _dim[ 8],_dim[ 9],_dim[10],_dim[11],
				 _dim[12],_dim[13],_dim[14],_dim[15]);

#pragma omp parallel default (shared) 
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			m[i] = (S)_M[i];
	}

    return m;

}

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

#endif // __MATRIX_H__
