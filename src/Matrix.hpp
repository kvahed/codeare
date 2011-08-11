/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
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
    COL, LIN, CHA, SET, ECO, PHS, REP, SEG, PAR, SLC, IDA, IDB, IDC, IDD, IDE, AVE, INVALID_DIM
};
#endif

#include "config.h"
#include "modules/OMP.hpp"

#include <complex>
#include <assert.h>

#include <iostream>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#ifdef HAVE_H5CPP_H
#include <H5Cpp.h>
using namespace H5;
#endif


#define ICE_SHRT_MAX 4095


/**
 * @brief raw data
 */
typedef std::complex<float> raw;


/**
 * Short test if the matrix is a vector.
 */
# define VECT(M) assert((M)->width() == 1 || (M)->height() == 1);

 
/**
 * Return absolute value.
 */
# define ABS(A) (A > 0 ? A : -A)


/**
 * Return minimum of two numbers
 */
# define MIN(A,B) (A > B ? A : B)


/**
 * Return maximum of two numbers 
 */
# define MAX(A,B) (A < B ? A : B)


/**
 * Return rounded value as in MATLAB.
 * i.e. ROUND( 0.1) ->  0
 *      ROUND(-0.5) -> -1
 *      ROUND( 0.5) ->  1
 */
# define ROUND(A) ( floor(A) + ((A - floor(A) >= 0.5) ? (A>0 ? 1 : 0) : 0))


/**
 * Some constants
 */
// PI
#ifndef PI
    # define PI  3.1415926535897931159979634685441851615906
#endif

// Gamma in Hz
#ifndef GAMMA
    #define GAMMA 4.2576e7
#endif

// Gamma in radians
#ifndef RGAMMA
    #define RGAMMA 267.513
#endif

#define KB 1024.0;
#define MB 1024.0 * 1024.0;
#define GB 1024.0 * 1024.0 * 1024.0;

/**
 * @brief   Matrix template.
 *          This class intends to offer a simple interface for handling
 *          MR data in a simple way. As of now it only support Siemens 
 *          access specifiers for direct input.
 *          The data is organised in a 16 member long array for dimensions
 *          and a template array for the data. The order is column-major.
 * 
 * @author  Kaveh Vahedipour
 * @date    Mar 2010
 */
template <class T>
class Matrix {
    

public:
    
    
    /**
     * @name Constructors and destructors
     *       Constructors and destructors
     */
    //@{
    
	
    /**
     * @brief           Contruct a 1^16 type matrix with (T)0.
     */
    inline              
    Matrix              ();
    
    
	/**
	 * @brief           Construct 16-dim matrix with dimension array
	 *
	 * @param  dim      All 16 Dimensions
	 */
	inline 
	Matrix              (const int* dim);
	
	
    /**
	 * @brief           Construct square 2D matrix
	 *
	 * @param  m        Rows & Columns
	 */
    inline              
    Matrix              (const int n) ;
    
    
    /**
	 * @brief           Construct 2D matrix
	 *
	 * @param  m        Rows
	 * @param  n        Columns
	 */
	inline 
	Matrix              (const int m, const int n);
	
    
    /**
	 * @brief           Construct 3D volume
	 *
	 * @param  m        Rows
	 * @param  n        Columns
	 * @param  k        Slices
	 */
	inline 
	Matrix              (const int m, const int n, const int k);
	

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
	inline 
	Matrix              (const int col, 
						 const int lin, 
						 const int cha,
						 const int set,
						 const int eco = 1,
						 const int phs = 1,
						 const int rep = 1,
						 const int seg = 1,
						 const int par = 1,
						 const int slc = 1,
						 const int ida = 1,
						 const int idb = 1,
						 const int idc = 1,
						 const int idd = 1,
						 const int ide = 1,
						 const int ave = 1);
	

    
    /**
     * @brief           Delete array containing data.
     */
    ~Matrix             ();


	/**
	 * @brief           Copy constructor
	 *
	 * @param  M        Right hand side
	 */
	inline 
	Matrix              (const Matrix<T> &M);

    //@}



    /**
     * @name            Import export functions for ICE access specifiers.
     *                  Ice access specifiers can be handled in one of the following ways.
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
    inline long         
    Import              (const IceAs ias);


    /**
     * @brief           Continue import from IceAs
     *                   
     * @param  ias      IceAs containing data
     * @param  pos      Import data starting at position pos of own repository
     * @return          Amount of data read
     */
    inline long         
    Import              (const IceAs ias, const long pos);


    /**
     * @brief           Export data to ias
     *                   
     * @param  ias      IceAs for data export
     * @return          Amount of data exported
     */
    inline long         
    Export              (IceAs* ias);


    /**
     * @brief           Partially export data to ias 
     * 
     * @param  ias      IceAs for data export
     * @param  pos      Export data starting at position pos of our repository
     */
    inline long         
    Export              (IceAs* ias, const long pos);
 
	#endif

    //@}

    
    /**
     * @name            Elementwise access
     *                  Elementwise access
     */
    
    //@{
    
    
    /**
     * @brief           Get pth element from data repository.
     *
     * @param  p        Requested position.
     * @return          Value at _M[p].
     */
    T                   
    operator[]          (const int p)                             const;
    
    
    /**
     * @brief           Reference to pth element from data repository.
     *
     * @param  p        Requested position.
     * @return          Reference to _M[p].
     */
    T                   
    &operator[]         (const int p)                              ;

    
    /**
     * @brief           Get pointer to data
     *  
     * @return          Data 
     */
    inline T*            
    Data                ()  const {
        return _M;
    }

    
    /**
     * @brief           Get element at position 
     *  
     * @param  pos      Position
     * @return          Value at _M[pos]
     */
    inline T            
    At                  (const int pos)  const {
        return _M[pos];
    }

    
    /**
     * @brief            Reference to value at position
     *  
     * @param  pos       Position
     * @return           Reference to _M[pos]
     */
    inline T&           
    At                  (const int pos) {
        return _M[pos];
    }

    
    /**
     * @brief           Get value in slice
     *  
     * @param  col      Column
     * @param  lin      Line
     * @return          Value at _M[col + _dim[LIN]*lin]
     */
    inline T            
    At                  (const int col, const int lin) const {

        return _M[col + _dim[COL]*lin ];

    }

    
    /**
     * @brief            Reference to value in slice
     *  
     * @param  col       Column
     * @param  lin       Line
     * @return           Reference to _M[col + _dim[LIN]*lin]
     */
    inline T&           
    At                  (int col, int lin) {

        return _M[col + _dim[COL]*lin ];

    }

    
    /**
     * @brief            Get value in volume
     *  
     * @param  col       Column
     * @param  lin       Line
     * @param  slc       Slice
     * @return           Value at _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc]
     */
    inline T            
    At                   (int col, int lin, int slc)  const {

        return _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc];

    }
    
    
    /**
     * @brief            Reference to value in volume
     *  
     * @param  col       Column
     * @param  lin       Line
     * @param  slc       Slice
     * @return           Reference to _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc]
     */
    inline T&            
    At                   (int col, int lin, int slc) {

        return _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc];

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
    At                   (const int col, 
						  const int lin, 
						  const int cha,
						  const int set,
						  const int eco,
						  const int phs = 0,
						  const int rep = 0,
						  const int seg = 0,
						  const int par = 0,
						  const int slc = 0,
						  const int ida = 0,
						  const int idb = 0,
						  const int idc = 0,
						  const int idd = 0,
						  const int ide = 0,
						  const int ave = 0) const {
		return _M [col+
				   lin*_dim[COL]+
				   cha*_dim[COL]*_dim[LIN]+
				   set*_dim[COL]*_dim[LIN]*_dim[CHA]+
				   eco*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]+
				   phs*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]+
				   rep*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]+
				   seg*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]+
				   par*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]+
				   slc*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]+
				   ida*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]+
				   idb*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]+
				   idc*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]+
				   idd*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]+
				   ide*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]*_dim[IDD]+
				   ave*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]*_dim[IDD]*_dim[IDE]];

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
    At                   (const int col, 
						  const int lin, 
						  const int cha,
						  const int set,
						  const int eco = 0,
						  const int phs = 0,
						  const int rep = 0,
						  const int seg = 0,
						  const int par = 0,
						  const int slc = 0,
						  const int ida = 0,
						  const int idb = 0,
						  const int idc = 0,
						  const int idd = 0,
						  const int ide = 0,
						  const int ave = 0) {
		return _M [col+
				   lin*_dim[COL]+
				   cha*_dim[COL]*_dim[LIN]+
				   set*_dim[COL]*_dim[LIN]*_dim[CHA]+
				   eco*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]+
				   phs*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]+
				   rep*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]+
				   seg*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]+
				   par*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]+
				   slc*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]+
				   ida*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]+
				   idb*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]+
				   idc*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]+
				   idd*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]+
				   ide*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]*_dim[IDD]+
				   ave*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]*_dim[IDD]*_dim[IDE]];

    }
    
	
	/**
	 * @brief           Create nxn identity matrix 
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn identity
	 */
	static Matrix<T> 
	Id                  (const int n);
	
	
	/**
	 * @brief           Create nxn matrix initialised with T(1.0)
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn ones
	 */
	static Matrix<T> 
	Ones                (const int n);
	
	
	/**
	 * @brief           Create nxn matrix initialised with T(1.0)
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn ones
	 */
	static Matrix<T> 
	Ones                (const int n, const int m);
	
	
	/**
	 * @brief           Create nxn matrix initialised with T(1.0)
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn ones
	 */
	static Matrix<T> 
	Ones                (const int n, const int m, const int l);
	
	
	/**
	 * @brief           Create nxn matrix initialised with T(0.0)
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn zeros
	 */
	static Matrix<T> 
	Zeros               (const int n);
	
	
	/**
	 * @brief           Create nxn matrix initialised with T(0.0)
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn zeros
	 */
	static Matrix<T> 
	Zeros               (const int n, const int m);
	
	
	/**
	 * @brief           Create nxn matrix initialised with T(0.0)
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn zeros
	 */
	static Matrix<T> 
	Zeros               (const int n, const int m, const int l);

	
	/**
	 * @brief           Reshape (MATLAB-like reshape). 
	 *                  All dims beyond lin are optional.
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
	Matrix <T>
	Reshape             (const int col, 
						 const int lin, 
						 const int cha = 1,
						 const int set = 1,
						 const int eco = 1,
						 const int phs = 1,
						 const int rep = 1,
						 const int seg = 1,
						 const int par = 1,
						 const int slc = 1,
						 const int ida = 1,
						 const int idb = 1,
						 const int idc = 1,
						 const int idd = 1,
						 const int ide = 1,
						 const int ave = 1) const {
		
		Matrix<T> res = (*this);
		res.reshape (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave);
		return res;
		
	}
	

	/**
	 * @brief           Reshape (MATLAB-like reshape). 
	 *                  All dims beyond lin are optional.
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
	void
	Reshape             (const int col, 
						 const int lin, 
						 const int cha = 1,
						 const int set = 1,
						 const int eco = 1,
						 const int phs = 1,
						 const int rep = 1,
						 const int seg = 1,
						 const int par = 1,
						 const int slc = 1,
						 const int ida = 1,
						 const int idb = 1,
						 const int idc = 1,
						 const int idd = 1,
						 const int ide = 1,
						 const int ave = 1) {
		
		int new_size = col * lin * cha * set * eco * phs * rep * 
			seg * par * slc * ida * idb * idc * idd * ide * ave;

		// Can't allow change of #elements
		assert (new_size == Size());

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

	}
	

	/**
	 * @brief           Get the element at position p of the vector, i.e. this(p).
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    T                  
    operator()          (const int p) const;

    
    /**
     * @brief           Get value of pth element of repository.
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    T&                 
    operator()          (const int p) ;

    
    /**
	 * @brief           Get value in slice
	 *
	 * @param  col      Column
	 * @param  lin      Line
	 * @return          Value at _M[col + _dim[LIN]*lin]
	 */
    T
    operator()          (const int col, const int lin) const {
        return _M[col + _dim[LIN]*lin ];
    }
    

    /**
	 * @brief           Reference to value in slice
	 *
	 * @param  col      Column
	 * @param  lin      Line
	 * @return          Reference to _M[col + _dim[LIN]*lin]
	 */
    T&                  
    operator()           (const int col, const int lin) {
        return _M[col + _dim[LIN]*lin ];
    }
    
    
    /**
     * @brief            Get value in volume
     *
     * @param  col       Column
     * @param  lin       Line
     * @param  slc       Slice
     * @return           Value at _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc]
     */
    T                  
    operator()           (const int col, const int lin, const int slc) const {
        return _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc];
    }
    
    
    /**
     * @brief            Reference to value in volume
     *
     * @param  col       Column
     * @param  lin       Line
     * @param  slc       Slice
     *
     * @return           Reference to _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc]
     */
    T&                 
    operator()           (const int col, const int lin, const int slc) {
           return _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc];
    }
    

	T&                 
    operator()           (const int col, 
						  const int lin, 
						  const int cha,
						  const int set,
						  const int eco = 0,
						  const int phs = 0,
						  const int rep = 0,
						  const int seg = 0,
						  const int par = 0,
						  const int slc = 0,
						  const int ida = 0,
						  const int idb = 0,
						  const int idc = 0,
						  const int idd = 0,
						  const int ide = 0,
						  const int ave = 0) { 

		return _M [col+
				   lin*_dim[COL]+
				   cha*_dim[COL]*_dim[LIN]+
				   set*_dim[COL]*_dim[LIN]*_dim[CHA]+
				   eco*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]+
				   phs*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]+
				   rep*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]+
				   seg*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]+
				   par*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]+
				   slc*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]+
				   ida*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]+
				   idb*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]+
				   idc*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]+
				   idd*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]+
				   ide*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]*_dim[IDD]+
				   ave*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]*_dim[IDD]*_dim[IDE]];
		
	}

	T
    operator()           (const int col, 
						  const int lin, 
						  const int cha,
						  const int set,
						  const int eco = 0,
						  const int phs = 0,
						  const int rep = 0,
						  const int seg = 0,
						  const int par = 0,
						  const int slc = 0,
						  const int ida = 0,
						  const int idb = 0,
						  const int idc = 0,
						  const int idd = 0,
						  const int ide = 0,
						  const int ave = 0) const { 
		
		return _M [col+
				   lin*_dim[COL]+
				   cha*_dim[COL]*_dim[LIN]+
				   set*_dim[COL]*_dim[LIN]*_dim[CHA]+
				   eco*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]+
				   phs*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]+
				   rep*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]+
				   seg*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]+
				   par*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]+
				   slc*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]+
				   ida*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]+
				   idb*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]+
				   idc*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]+
				   idd*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]+
				   ide*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]*_dim[IDD]+
				   ave*_dim[COL]*_dim[LIN]*_dim[CHA]*_dim[SET]*_dim[ECO]*_dim[PHS]*_dim[REP]*_dim[SEG]*_dim[PAR]*_dim[SLC]*_dim[IDA]*_dim[IDB]*_dim[IDC]*_dim[IDD]*_dim[IDE]];

	}

    //@}
    


    /**
     * @name            Partial copy functions.
     *                  Functions for access to parts of data.
     */
    
    //@{
    

    /**
     * @brief            Operates only on inner 2D: Get a Row of data
     *  
     * @param  r         Row
     * @return           Copy data into new vector
     */
    Matrix <T>           
    Row                  (int r) const;
	
    
    /**
     * @brief            Operates only on inner 2D: Get a Row of data
     *  
     * @param  c         Column
     * @return           Copy data into new vector
     */
    Matrix <T>           
    Column               (int c) const;
	
    /**
     * @brief            Operates only on inner 3D: Get a slice of data
     *  
     * @param  s         Slice
     * @return           Copy data into new matrix and return
     */
    Matrix <T>           
    Slice                (int s) const;
	
    /**
     * @brief            Operates only on inner 3D: Get a slice of data
     *  
     * @param  s         
     * @return           Copy data into new matrix and return
     */
    Matrix <T>           
    Volume               (int v) const;
	

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
    inline int                 
    height              () const {return _dim[0];}
    
    
    /**
     * @brief           Get refernce to number of rows, i.e. tmp = size(this); tmp(1).
     *
     * @return          Number of rows.
     */
    inline int&
    height              () {return _dim[0];}
    
    
    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    inline int                 
    width               () const {return _dim[1];}
    
    
    /**
     * @brief           Get reference number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    inline int&
    width               ()  {return _dim[1];}
    
    
    /**
     * @brief           Get number of rows, i.e. tmp = size(this); tmp(1).
     *
     * @return          Number of rows.
     */
    int                 
    m                   () const {return _dim[0];}


    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    int                 
    n                   () const {return _dim[1];}


    /**
     * @brief           Get size a given dimension.
     *
     * @return          Number of rows.
     */
    inline int          
    Dim                 (const int i)                                const {return _dim[i];}
    
    
    /**
     * @brief           Get reference to size a given dimension.
     *
     * @return          Number of rows.
     */
    inline int&          
    Dim                 (const int i)                                 {return _dim[i];}
    
    
    /**
     * @brief           Get dimension array
     *
     * @return          All dimensions
     */
    inline const int*   
    Dim                 ()                                     const {return _dim;}
    

    /**
     * @brief           Reset all dimensions to values in dim
	 *
	 * @param  dim      New dimensions
     */
    inline void         
    Dim                 (const int* dim)                                     const {
        for (int i = 0; i<INVALID_DIM; i++)
            _dim[i] = dim[i];
    }
    
    
    /**
     * @brief           Resize to dims, reallocate and zero repository. 
	 *                  Needs to becalled after any resize operation on dimensios. 
	 *                  (i.e. M.Dim(2) = 10; Reset();)
     *
     * @param  dim      New dimensions
     */
    inline void         
    Reset               (const int* dim)                                      {

    	for (int i = 0; i < INVALID_DIM; i++)
            _dim[i] = dim[i];

        if (nb_alloc) {
            free (_M);
            nb_alloc--;
        }

        _M = (T*) malloc (Size()*sizeof(T));
        nb_alloc++;

		Zero();

    }
    

    /**
     * @brief           Purge data and free RAM.
     */
    inline void         
    Clear               ()                                      {

    	for (int i = 0; i < INVALID_DIM; i++)
            _dim[i] = 1;

        if (nb_alloc) {
            free (_M);
            nb_alloc--;
        }

    }
    

    /**
     * @brief           Reset. i.e. reallocate and set all fields = T(0)
     */
    inline void         
    Reset               ()                                      {

        if (nb_alloc) {
            free (_M);
            nb_alloc--;
        }
		
		_M = (T*) malloc (Size() * sizeof (T));
        nb_alloc++;

		Zero();

    }
    

    /**
     * @brief           Reset. i.e. Set all fields = T(0)
     */
    inline void         
    Zero               ()                                      {

		for (int i = 0; i < Size(); i++)
			_M[i] = (T) 0;

    }
    

    /**
     * @brief           Get the number of matrix cells, i.e. dim_0*...*dim_16.
     *
     * @return          Number of cells.
     */
    long                
    Size                ()                                    const;
    
    
    /**
     * @brief           Check if we are XD (i.e. X dimensions > 1)
     *
	 * @param  d        Dimensions
     * @return          XD matrix?
     */
    bool                
    IsXD                (const int d)                         const;
    
    
    /**
     * @brief           Check if we are 2D (i.e. COL, LIN)
     *
     * @return          2D matrix?
     */
    bool                
    Is1D                ()                                    const;
    
    
    /**
     * @brief           Check if we are 2D (i.e. COL, LIN)
     *
     * @return          2D matrix?
     */
    bool                
    Is2D                ()                                    const;
    
    
    /**
     * @brief           Check if we are 3D (i.e. COL, LIN, SLC)
     *
     * @return          3D matrix?
     */
    bool                
    Is3D                ()                                    const;
    
    
    /**
     * @brief           Check if we are 4D (i.e. COL, LIN, SLC)
     *
     * @return          3D matrix?
     */
    bool                
    Is4D                ()                                    const;
    
    
    /**
     * @brief           Get the number of matrix cells, i.e. Size * sizeof(T).
     *
     * @return          Size in RAM in bytes.
     */
    int                 
    SizeInRAM           ()                                    const;

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
    Matrix<T>           
    operator=           (Matrix<T> &M);
    
    
    /**
     * @brief           Assignment operator. i.e. this = m.
     *
     * @param  M        The assigned matrix.
     */
    Matrix<T>           
    operator=           (const Matrix<T> &M);
    
    
    /**
     * @brief           Matrix product. i.e. this * M.
     *
     * @param  M        The factor.
     */
    Matrix<T>           
    operator->*         (Matrix<T> &M);
    
    
    /**
     * @brief           Elementwise substruction of two matrices
     *
     * @param  M        Matrix substruent.
     */
    Matrix<T>           
    operator-           (Matrix<T> &M);
    
    
    /**
     * @brief           Elementwise substruction all elements by a scalar
     *
     * @param  s        Scalar substruent.
     */
    Matrix<T>           
    operator-           (T s);
    
    
    /**
     * @brief           Negate or substruct from 0. i.e. 0 - this.
	 *
	 * @return          Negation
     */
    Matrix<T>           
    operator-           ();
    
    
    /**
     * @brief           Elementwise addition of two matrices
     *
     * @param  M        Matrix substruent.
     */
    Matrix<T>           
    operator+           (Matrix<T> &M);
    
    
    /**
     * @brief           Elementwise addition iof all elements with a scalar
     *
     * @param  s        Scalar substruent.
     */
    Matrix<T>           
    operator+           (T s);
    
    
    /**
     * @brief           Add to 0
	 *
	 * @return          M+0;
     */
    Matrix<T>           
    operator+           ();
    
    
    /**
     * @brief           Transposition / Complex conjugation. i.e. this'.
	 *
	 * @return          Matrix::tr()
     */
    Matrix<T>           
    operator!           () const;
    
    
    /**
     * @brief           Return a matrix with result[i] = (m[i] ? this[i] : 0).
	 *
	 * @param           The operand
	 * @return          Cross-section or zero
     */
    Matrix<T>           
    operator&           (Matrix<bool> &M);
    
    
    /**
     * @brief           Scalar equality. result[i] = (this[i] == m).
     *
     * @param  s        Comparing scalar.
     * @return          Matrix of true where elements are equal s and false else.
     */
    Matrix<bool>        
    operator==          (T s);
    
    
    /**
     * @brief           Scalar inequality. result[i] = (this[i] != m). i.e. this ~= m
     *
     * @param  s        Comparing scalar.
	 * @param           Matrix of false where elements are equal s and true else.
     */
    Matrix<bool>        
    operator!=          (T s);
    
    
    /**
     * @brief           Scalar greater comaprison, result[i] = (this[i] > m). i.e. this > m
     *
     * @param  s        Comparing scalar.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator>           (T s);
    
    
    /**
     * @brief           Scalar greater or equal comparison. result[i] = (this[i] >= m). i.e. this >= m
     *
     * @param  s        Comparing scalar.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator>=          (T s);
    
    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] <= m). i.e. this <= m
     *
     * @param  s        Comparing scalar.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator<=          (T s);
    
    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] < m). i.e. this < m
     *
     * @param  s        Comparing scalar.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator<           (T s);
    
    
    /**
     * @brief           Elementwise equality, result[i] = (this[i] == m[i]). i.e. this == m
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator==          (Matrix<T> M);
    
    
    /**
     * @brief           Elementwise equality, result[i] = (this[i] != m[i]). i.e. this ~= m
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator!=          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] >= m[i]). i.e. this >= m
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator>=          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] <= m[i]). i.e. this <= m.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator<=          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] > m[i]). i.e. this > m.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator>           (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] < m[i]). i.e. this < m.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator<           (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (m[i] || this[i] ? 1 : 0). i.e. this | m.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<T>           
    operator||          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (m[i] && this[i] ? 1 : 0). i.e. this & m.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<T>           
    operator&&          (Matrix<T> &M);


    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<T>           
    operator^           (float p);
    
	//friend Matrix<T> operator+(Matrix<T> &a, Matrix<T> &b){return a+b;};
    //@}





    /**
     * @name Data manipulation functions.
     *       Data manipulation functions.
     */
    
    //@{
    
    /**
     * @brief           Fill with random data
     */
    inline void         
	Random               ();    


    /**
     * @brief           Elementwise multiplication. i.e. this .* M.
     *
     * @param  M        Factor matrix.
	 * @retrun          Result
     */
    Matrix<T>           
    operator*           (Matrix<T> &M);


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
	 * @retrun          Result
     */
    Matrix<T>           
    operator*=           (Matrix<T> &M);
    
    
    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Factor scalar.
	 * @retrun          Result
     */
    Matrix<T>           
    operator*=           (T s);
    
    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
	 * @retrun          Result
     */
    Matrix<T>           
    operator+=           (Matrix<T> &M);
    
    
    /**
     * @brief           ELementwise addition with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
	 * @retrun          Result
     */
    Matrix<T>           
    operator+=           (T s);
    
    
    /**
     * @brief           Elementwise multiplication with a scalar. i.e. this * m.
	 *
	 * @param  s        Factor scalar
	 * @return          Result
     */
    Matrix<T>           
    operator*           (T s);

    
    /**
     * @brief           Elelemtwise division by M.
     *
     * @param  M        The divisor.
	 * @return          Result
     */
    Matrix<T>           
    operator/           (Matrix<T> &M);

    
    /**
     * @brief           Elementwise division by scalar. i.e. this * m.
	 *
     * @param  s        The divisor.
	 * @return          Result
     */
    Matrix<T>           
    operator/           (T s);
    
    //@}
    
    /**
     * @name Display and IO functions.
     *       Display and IO functions.
     */
    
    //@{
    

    /**
     * @brief           Print contents to output stream.
     *
     * @param  os       The output stream.
     * @return          The output stream.
     */
    std::ostream&       
    Print               (std::ostream &os) const;
    

    /**
     * @brief           Print contents to output stream.
     *
     * @param  os       The output stream.
     * @return          The output stream.
     */
    std::ostream&
    operator<<          (std::ostream &os);


    /**
     * @brief           Primitive dump column-major to file.
     * 
     * @param  fname    File name.
     * @return          Success.
     */
    bool                
    pdump               (std::string fname);
    

    /**
     * @brief           Dump to HDF5 file.
     * 
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    h5dump                (std::string fname, std::string dname = "", std::string dloc = "/");
    

    /**
     * @brief           Read from HDF5 file.
     *
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    h5read              (std::string fname, std::string dname = "", std::string dloc = "/");
    
    /**
     * @brief           Dump to MATLAB file.
     * 
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    mxdump              (std::string fname, std::string dname = "", std::string dloc = "/");
    

    /**
     * @brief           Read from MATLAB file.
     *
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    mxread              (std::string fname, std::string dname = "", std::string dloc = "/");
    
    /**
     * @brief           Dump to HDF5 file.
     * 
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    dump                (std::string fname, std::string dname = "", std::string dloc = "/", std::string fmt = "h5");
    

    /**
     * @brief           Read from HDF5 file.
     *
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    read                (std::string fname, std::string dname = "", std::string dloc = "/", std::string fmt = "h5");
    
	/**
	 * @brief 
	 */
	bool
	rawread             (const std::string fname, const std::string version);

	/**
	 * @brief 
	 */
	bool
	rsadjust            (const std::string fname);

	/**
	 * @brief           Skip lines in file
	 *
	 * @param  file     Position pointer in stream
	 * @param  lines    Lines to skip
	 * @return          New position
	 */
	std::fstream&
	SkipLines           (const std::ifstream& file, const unsigned int num);

    //@}
    
    
    
    /**
     * @name            Other functions.
     *                  Other functions.
     */
    
    //@{
    
    
    /**
     * @brief           Greatest element of the matrix. i.e. max(M0).
     *
     * @return          The scalar value of the greatest element.
     */
    T                   
    Max                 ();    
    
    
    /**
     * @brief           Smallest element of the matrix. i.e. min(M0).
     *
     * @return          The scalar value of the smallest element.
     */
    T                   
    Min                 ();
    
    
    /**
     * @brief           Get maximum absolute value in the matrix
     *
     * @return          Maximum value
     */
    T                   
    Maxabs              ();
    
    
    /**
     * @brief           Get minimum absolute value in the matrix
     *
     * @return          Maximum value
     */
    T                   
    Minabs              ();
    
    /**
     * @brief           Matrix Product.
     *
     * @return          Product of this and M.
     */
    Matrix<T>           
    prod                (Matrix<T> &M);
    
    /**
     * @brief           Matrix Product with transpose / complex conjugate transpose.
     *
     * @return          Product of this and M.
     */
    Matrix<T>           
    prodt               (Matrix<T> &M);
    
    /**
     * @brief           Transposition.
     *
     * @return          The transposed matrix
     */
    Matrix<T>           
    tr   ()             const;
    
    /**
     * @brief           Fourier transform over all dimensions.
     *
     * @return          Fourier transform.
     */
    Matrix<T>           
    FFT                 ()            const;
    
    /**
     * @brief           Inverse Fourier transform over all dimensions
     *
     * @return          Inverse Fourier transform.
     */
    Matrix<T>           
    IFFT                ()            const;
    
    /**
     * @brief           MATLAB-like fftshift; 
	 *                  i.e. shift zero component to centre of volume for viewing
     *
     * @return          FFT shift.
     */
    Matrix<T>           
    FFTShift            (const int d = 0)            const;
    
    /**
     * @brief           MATLAB-like fftshift; 
	 *                  i.e. shift zero component to centre of volume for viewing
     *
     * @return          FFT shift.
     */
    Matrix<T>           
    IFFTShift           (const int d = 0)            const;
    

    /**
     * @brief           ND Hann window
     *
     * @return          FFT shift.
     */
    Matrix<T>           
    HannWindow          (const int d = 0)            const;
    

    /**
     * @brief           Sum of squares. 
     *
	 * @param  d        Dimension to eliminate. If not given, it is done over the outermost
     * @return          FFT shift.
     */
	Matrix<T>
    SOS                 (const int d = 0)           const;
    

    /**
     * @brief           Squeeze dimensions
     */
	void
    Squeeze             ();
    

    /**
     * @brief           Squeeze dimensions
     */
	Matrix<T>
    Squeeze             () const ;
    

	/**
	 * @brief           Highest occupied dimension
	 *
	 * @return 
	 */
	unsigned short
	HDim                () const;


	/**
	 * @brief           Print dimension to std out
	 *
	 * @retutn 
	 */
	void
	PrintDims           () const;


    /**
     * @brief           General inversion.
     *
     * @return          success
     */
	int 
    Inv   ()            const;
    
    /**
     * @brief           Moore-Penrose pseudo-inverse.
     *
     * @return          Matrix
     */
	Matrix<T> 
    Pinv   ();
    
    /**
     * @brief           Euclidean norm.
     *
     * @return          Norm
     */
	T
    Norm ()             const;
    
    /**
     * @brief           Dot product, complex, conjugate first vector
     *
     * @param  M        Factor      
     */
	T
    dotc (Matrix<T>& M) const;
    
    //@}


    /**
     * @name Lapack functions.
     *       Only available when lapack detected.
     */
    
    //@{
    

    /**
     * @brief           Compute eigen values with Lapack
     *
	 * @param  ev       Vector containing the computed eigenvalues
	 * @param  cv       Compute also left hand eigen vectors
	 * @param  lev      Left hand eigen vectors 
	 * @param  rev      right hand eigen vectors 
     * @return          Feedback from Lapack operation
     */
    inline int
	EIG                 (const bool cv, Matrix<raw>* ev, Matrix<T>* lev, Matrix<T>* rev);
    

    /**
     * @brief           Compute singular value decomposition with lapack.
     *
	 * @param  cm       Compute left and right matrices U and V**T
	 * @param  lsv      Left hand singular vectors.
	 * @param  rsv      Right hand singular vectors.
	 * @param  sv       Sorted singular values.
     * @return          Info from Lapack operation.
     */
    inline int
	SVD                 (const bool cm, Matrix<T>* lsv, Matrix<T>* rsv, Matrix<T>* sv);
    
	
    //@}
    
    

private:
    
    int                 _dim[INVALID_DIM]; /// Dimnesions
    T*                  _M;                /// Data repository
	int                 nb_alloc;

    /**
     * @brief           Matrix Product with BLAS.
     *
	 * @param  M        Multiplie with
	 * @param  transb   Transpose M
	 *
     * @return          Product of this and M.
     */
    Matrix<T>           
    GEMM                (Matrix<T> &M, char transb);
    
};


template <class T> Matrix<T> 
Matrix<T>::Volume (int s) const {
    
	assert (Is4D());
    
    Matrix<T> res;

	for (int j = 0; j < 3; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	int nc = _dim[0]*_dim[1]*_dim[2];

	memcpy (&res[0], &_M[s * nc], nc * sizeof(T));

	return res;

}


template <class T> Matrix<T> 
Matrix<T>::Slice (int s) const {
    
	assert (Is3D());
    
    Matrix<T> res;

	for (int j = 0; j < 2; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	int nc = _dim[0]*_dim[1];

	memcpy (&res[0], &_M[s * nc], nc*sizeof(T));

	return res;

}


template <class T> Matrix<T> 
Matrix<T>::Row (int r)  const {

	assert (Is2D());
    
    Matrix<T> res;

	res.Dim(0) = _dim[1];
	res.Reset();

	for (int i = 0; i < _dim[1]; i++)
		res[i] = _M[r + i*_dim[0]];

	return res;

}


template <class T> Matrix<T> 
Matrix<T>::Column (int c) const {
    
    Matrix<T> res;

	res.Dim(0) = _dim[0];
	res.Reset();

	memcpy (&res[0], _M[c*_dim[0]], _dim[0] * sizeof(T));

	return res;

}


template <class T> long 
Matrix<T>::Size() const {
    
    long size = 1;
    
    for (int i = 0; i < INVALID_DIM; i++)
        size *= _dim[i];
    
    return size;
    
}


template <class T> inline int 
Matrix<T>::SizeInRAM() const {
    
    return Size() * sizeof(T);
    
}


template <> inline short 
Matrix<short>::Max() {
	
    short max = _M[0];
	
    for (int i = 0; i < Size(); i++)
        if (_M[i] > max)
            max = _M[i];
	
    return max;
	
}


template <> inline raw
Matrix<raw>::Max() {
		
		raw   max = raw(0.0,0.0);
		float tmp =  0.0;
		
		for (int i = 0; i < Size(); i++) {
				float abs = sqrt(_M[0].real()*_M[0].real() + _M[0].imag()*_M[0].imag());
				if (abs > tmp) {
						tmp = abs;
						max = _M[i];
				}
		}
		return max;
	
}


template <class T>
T  Matrix<T>::Maxabs() {

    T max = fabs(_M[0]);

    for (int i = 0; i < Size(); i++)
        if (fabs(_M[i]) > max)
            max = fabs(_M[i]);

    return max;

}


template <class T>
T Matrix<T>::Min() {

    T min = _M[0];

    for (int i = 0; i < Size(); i++)
        if (_M[i] < min)
            min = _M[i];

    return min;

}


template <class T>
T  Matrix<T>::Minabs() {

    T old = fabs(_M[0]);

    for (int i = 0; i < Size(); i++)
        if (fabs(_M[i]) < old)
            old = fabs(_M[i]);

    return old;

}

template <class T>
Matrix<T> Matrix<T>::tr() const {

    Matrix<T> res (_dim);
	
	long tmp   = res.Dim(0);
	res.Dim(0) = res.Dim(1);
	res.Dim(1) = tmp;
	
    for (int i = 0; i < res.Dim(0); i++)
        for (int j = 0; j < res.Dim(1); j++)
			if (typeid (T) == typeid (double))
				res.At(i,j) =      At(j,i);  // Transpose
			else
				res.At(i,j) = conj(At(j,i)); // Conjugate transpose

    return res;

}

template<>
inline void Matrix<raw>::Random () {
	
	srand (time(NULL));

	for (int i = 0; i < Size(); i++)
		_M[i] = raw ((float) rand() / (float) RAND_MAX*2-1, (float) rand() / (float) RAND_MAX*2-1);
	
}
    
template<>
inline void Matrix<double>::Random () {

	srand (time(NULL));

	for (int i = 0; i < Size(); i++)
		_M[i] = (double) rand() / (double) RAND_MAX*2-1;

}
    
template<>
inline void Matrix<short>::Random () {

	srand (time(NULL));

	for (int i = 0; i < Size(); i++)
		_M[i] = (short) 12 * (double)rand() / (double)RAND_MAX*2-1;

}

    
template <class T> 
inline bool Matrix<T>::IsXD (const int d) const {

	int l = 0;

	for (int i = 0; i < INVALID_DIM; i++)
		if (_dim[i] > 1) l++;

	return (l == d);

}

template <class T> 
inline bool Matrix<T>::Is1D () const {
	
	return IsXD(1);

}

template <class T> 
inline bool Matrix<T>::Is2D () const {
	
	return IsXD(2);

}

template <class T> 
inline bool Matrix<T>::Is3D () const {
	
	return IsXD(3);

}

template <class T> 
inline bool Matrix<T>::Is4D () const {
	
	return IsXD(4);

}

#include "Matrix_Constructors.hpp"
#include "Matrix_IO.cpp"
#include "Matrix_Lapack.cpp"
#include "Matrix_Operators.cpp"
#include "Matrix_BLAS.cpp"
#include "Matrix_ICE.cpp"
#include "Matrix_FFT.hpp"

#endif // __MATRIX_H__
