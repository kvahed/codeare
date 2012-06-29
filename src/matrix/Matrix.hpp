/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
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
#include "modules/OMP.hpp"
#include "Complex.hpp"

#include "cycle.h"            // FFTW cycle implementation

#include <assert.h>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <vector>
#include <ostream>

#define ICE_SHRT_MAX 4095


#ifdef HAVE_MAT_H
#include "mat.h"
#endif

#include "Ptr.h"


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
 * Some constants
 */
// PI
#ifndef PI
    # define PI 3.141592653589793238462643383279502884197169399375105820974944592
#endif

#ifndef TWOPI
    #define TWOPI 6.283185307179586476925286766559005768394338798750211641949889185
#endif

// Gamma in Hz
#ifndef GAMMA
    #define GAMMA 42.576
#endif

// Gamma in radians
#ifndef RGAMMA
    #define RGAMMA 267.513
#endif

#define KB 1024.0;
#define MB 1024.0 * 1024.0;
#define GB 1024.0 * 1024.0 * 1024.0;
#define TB 1024.0 * 1024.0 * 1024.0 * 1024.0;

/**
 * @brief   Matrix template.<br/>
 *          This class intends to offer a simple interface for handling
 *          MR data in a simple way. <br/>As of now it only support Siemens 
 *          access specifiers for direct input.<br/>
 *          The data is organised in a 16 member long array for dimensions
 *          and a template array for the data. <br/>The order is column-major.
 * 
 * @author  Kaveh Vahedipour
 * @date    Mar 2010
 */
template <class T>
class Matrix  : public SmartObject {
    

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
	Matrix              (const size_t* dim);
	
	
	/**
	 * @brief           Construct 16-dim matrix with dimension and resolution arrays
	 *
	 * @param  dim      All 16 Dimensions
	 * @param  res      All 16 Resolutions
	 */
	inline 
	Matrix              (const size_t* dim, const float* res);
	
	
    /**
	 * @brief           Construct square 2D matrix
	 *
	 * @param  n        Rows & Columns
	 */
    inline              
    Matrix              (const size_t& n) ;
    
    
    /**
	 * @brief           Construct 2D matrix
	 *
	 * @param  m        Rows
	 * @param  n        Columns
	 */
	inline 
	Matrix              (const size_t& m, const size_t& n);
	
    
    /**
	 * @brief           Construct 3D volume
	 *
	 * @param  m        Rows
	 * @param  n        Columns
	 * @param  k        Slices
	 */
	inline 
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
	inline 
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
	 * @param  M        Right hand side
	 */
	inline 
	Matrix              (const Matrix<T> &M);


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
     * @brief           Get pth element from data repository.
     *
     * @param  p        Requested position.
     * @return          Value at _M[p].
     */
    T                   
    operator[]          (const size_t& p) const;
    
    
    /**
     * @brief           Reference to pth element from data repository.
     *
     * @param  p        Requested position.
     * @return          Reference to _M[p].
     */
    T                   
    &operator[]         (const size_t& p);

    
    /**
     * @brief           Get pointer to data
     *  
     * @return          Data 
     */
    inline const T*            
    Data                (const size_t pos = 0)  const {
        return &(_M.at(pos));
	}

    
    /**
     * @brief           Get pointer to data
     *  
     * @return          Data 
     */
    inline std::vector<T>&            
    Dat                 ()  {
        return _M;
	}

    
    /**
     * @brief           Get pointer to data
     *  
     * @return          Data 
     */
    inline std::vector<T>            
    Dat                 ()  const {
        return _M;
	}

    
    /**
     * @brief           Get element at position 
     *  
     * @param  pos      Position
     * @return          Value at _M[pos]
     */
    inline T            
    At                  (const size_t& pos) const {

        return _M[pos];

    }

    

    /**
     * @brief            Reference to value at position
     *  
     * @param  pos       Position
     * @return           Reference to _M[pos]
     */
    inline T&           
    At                  (const size_t& pos) {

        return _M[pos];

    }


    
    /**
     * @brief           Get value in slice
     *  
     * @param  col      Column
     * @param  lin      Line
     * @return          Value at _M[col + _dim[COL]*lin]
     */
    inline T            
    At                  (const size_t& col, const size_t& lin) const {

        return _M[col + _dim[COL]*lin ];

    }

    

    /**
     * @brief            Reference to value in slice
     *  
     * @param  col       Column
     * @param  lin       Line
     * @return           Reference to _M[col + _dim[COL]*lin]
     */
    inline T&           
    At                  (const size_t& col, const size_t& lin) {

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
    At                   (const size_t& col, const size_t& lin, const size_t& slc)  const {

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
    At                   (const size_t& col, const size_t& lin, const size_t& slc) {

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
	 * @brief           Reshape (MATLAB-like reshape).<br/>
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
	inline Matrix <T>
	Reshape             (const size_t& col, 
						 const size_t& lin, 
						 const size_t& cha = 1,
						 const size_t& set = 1,
						 const size_t& eco = 1,
						 const size_t& phs = 1,
						 const size_t& rep = 1,
						 const size_t& seg = 1,
						 const size_t& par = 1,
						 const size_t& slc = 1,
						 const size_t& ida = 1,
						 const size_t& idb = 1,
						 const size_t& idc = 1,
						 const size_t& idd = 1,
						 const size_t& ide = 1,
						 const size_t& ave = 1) const {
		
		Matrix<T> res = (*this);
		res.reshape (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave);
		return res;
		
	}
	

	/**
	 * @brief           Reshape (MATLAB-like reshape).<br/> 
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
	inline void
	Reshape             (const size_t& col, 
						 const size_t& lin, 
						 const size_t& cha = 1,
						 const size_t& set = 1,
						 const size_t& eco = 1,
						 const size_t& phs = 1,
						 const size_t& rep = 1,
						 const size_t& seg = 1,
						 const size_t& par = 1,
						 const size_t& slc = 1,
						 const size_t& ida = 1,
						 const size_t& idb = 1,
						 const size_t& idc = 1,
						 const size_t& idd = 1,
						 const size_t& ide = 1,
						 const size_t& ave = 1) {
		
		size_t new_size = col * lin * cha * set * eco * phs * rep * 
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
	 * @brief          Cast operator
	 *
	 * @return         Cast if possible
	 */
	template<class S> operator Matrix<S> () const;


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
	 * @param  col      Column
	 * @param  lin      Line
	 * @return          Value at _M[col + _dim[COL]*lin]
	 */
    T
    inline 
	operator()          (const size_t& col, const size_t& lin) const {

        return this->At(col, lin);

    }
    

    /**
	 * @brief           Reference to value in slice
	 *
	 * @param  col      Column
	 * @param  lin      Line
	 * @return          Reference to _M[col + _dim[COL]*lin]
	 */
    inline T&                  
    operator()           (const size_t& col, const size_t& lin) {

        return this->At(col, lin);

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
    operator()           (const size_t& col, const size_t& lin, const size_t& slc) const {

        return this->At(col, lin, slc);

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
    inline T&                 
    operator()           (const size_t& col, const size_t& lin, const size_t& slc) {

		return this->At(col, lin, slc);

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
     *                  Friend oprators
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
	inline friend Matrix<T>    
	operator*  (const double& s, const Matrix<T>& m) { 
		return   m * s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator*  (const float& s, const Matrix<T> &m) { 
		return   m * s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator*  (const short& s, const Matrix<T> &m) { 
		return   m * s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator*  (const long& s, const Matrix<T> &m) { 
		return   m * s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator*  (const cxfl& s, const Matrix<T> &m) { 
		return   m * s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator*  (const cxdb& s, const Matrix<T> &m) { 
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
	inline friend Matrix<T>    
	operator+  (const double& s, const Matrix<T> &m) { 
		return   m + s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator+  (const float& s, const Matrix<T> &m) { 
		return   m + s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator+  (const short& s, const Matrix<T> &m) { 
		return   m + s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator+  (const long& s, const Matrix<T> &m) { 
		return   m + s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator+  (const cxfl& s, const Matrix<T> &m) { 
		return   m + s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator+  (const cxdb& s, const Matrix<T> &m) { 
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
	inline friend Matrix<T>    
	operator-  (const double& s, const Matrix<T> &m) { 
		return   m - s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator-  (const float& s, const Matrix<T> &m) { 
		return   m - s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator-  (const short& s, const Matrix<T> &m) { 
		return   m - s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator-  (const long& s, const Matrix<T> &m) { 
		return   m - s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator-  (const cxfl& s, const Matrix<T> &m) { 
		return   m - s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator-  (const cxdb& s, const Matrix<T> &m) { 
		return   m - s;
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator/  (const double& s, const Matrix<T> &m) { 

		Matrix<T> res = m;
		size_t i = res.Size();

		while (i--)
			res[i] = (res[i] != 0) ? s/res[i] : 0.0;
		
		return res;
		
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	inline friend Matrix<T>    
	operator/  (const float& s, const Matrix<T> &m) { 

		Matrix<T> res = m;
		size_t i = res.Size(); 

		while (i--)
			res[i] = (res[i] != 0) ? s/res[i] : 0.0;

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
	operator== (const T& s, const Matrix<T>& m) {
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
	operator>= (const T& s, const Matrix<T>& m) {
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
	operator<= (const T& s, const Matrix<T>& m) {
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
	operator!= (const T& s, const Matrix<T>& m) {
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
	operator>  (const T& s, const Matrix<T>& m) {
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
	operator<  (const T& s, const Matrix<T>& m) {
		return   m >  s;
	}


	/**
	 * @brief           Elementwise equality with scalar (lhs)
	 *
	 * @param  mb       Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          T+M
	 */
	inline friend Matrix<T>    
	operator&  (const Matrix<bool>& mb, const Matrix<T>& m) {
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
    Height              () const {return _dim[0];}
    
    
    /**
     * @brief           Get refernce to number of rows, i.e. tmp = size(this); tmp(1).
     *
     * @return          Number of rows.
     */
    inline size_t&
    Height              () {return _dim[0];}
    
    
    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    inline size_t                 
    Width               () const {return _dim[1];}
    
    
    /**
     * @brief           Get reference number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    inline size_t&
    Width               ()  {return _dim[1];}
    
    
    /**
     * @brief           Get resolution a given dimension.
     *
	 * @param   i       Dimension
     * @return          Resolution .
     */
    inline float          
    Res                 (const size_t& i)                                const {return _res[i];}
    
    
    /**
     * @brief           Rresolution a given dimension.
     *
	 * @param   i       Dimension
     * @return          Resolution
     */
    inline float&          
    Res                 (const size_t& i)                                 {return _res[i];}
    


	/**
	 * @brief           Resolution array
	 *
	 * @return          All resolutions
	 */
	const float*
	Res                 () const {
		return &_res[0];
	};
	

    
    /**
     * @brief           Get size a given dimension.
     *
	 * @param   i       Dimension
     * @return          Dimension
     */
    inline size_t          
    Dim                 (const size_t& i)                                const {return _dim[i];}
    
    
    /**
     * @brief           Get reference to size a given dimension.
     *
     * @return          Number of rows.
     */
    inline size_t&          
    Dim                 (const size_t& i)                                 {return _dim[i];}
    
    
    /**
     * @brief           Get dimension array
     *
     * @return          All dimensions
     */
    inline const size_t*   
    Dim                 ()                                     const {return _dim;}
    

    /**
     * @brief           Get size a given dimension.
     *
     * @return          Number of rows.
     */
    inline size_t          
    Dim                 (const int& i)                                const {return _dim[i];}
    
    
    /**
     * @brief           Get reference to size a given dimension.
     *
     * @return          Number of rows.
     */
    inline size_t&          
    Dim                 (const int& i)                                 {return _dim[i];}
    
    
    /**
     * @brief           Reset all dimensions to values in dim
	 *
	 * @param  dim      New dimensions
     */
    inline void         
    Dim                 (const size_t* dim)    {

		for (size_t i = 0; i < INVALID_DIM; i++)
			_dim [i] = dim [i];

		Reset();

    }
    

	/**
	 * @brief           Expand matrix by increasing highest dimension (Concatenation)
	 *
	 * @param  dim      Dimension to be expanded
	 * @param  n        Expand by n x current size (default 1)
	 */
	inline void
	Expand              (const size_t dim, const size_t& n = 1)                      {

		_dim[dim] += n;
		_M.resize(Size());
		
	}

    
	/**
	 * @brief           Expand matrix by increasing highest dimension (Concatenation)
	 *
	 * @param  dim      Dimension to be expanded
	 * @param  n        Expand by n x current size (default 1)
	 */
	inline void
	PushBack            (const Matrix<T> M)                    {

		size_t nd = 0;
		size_t os = Size(); // old size 
		
		// highest non-singlton dimension
		for (size_t i = 0; i < INVALID_DIM; i++)
			nd  = (M.Dim(i) > 1) ? i : nd;

	   
		
		_dim[nd] += 1;
		
		_M.resize(Size());

		
		
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

		_M.resize(Size(), T(0));

    }
    

    /**
     * @brief           Resize
     */
    inline void         
    Resize              () {

		_M.resize(Size(), T(0));

    }
    

    /**
     * @brief           Resize to dims, reallocate and zero repository. 
	 *                  Needs to becalled after any resize operation on dimensios. 
	 *                  (i.e. M.Dim(2) = 10; Reset();)
     *
     * @param  dim      New dimensions
     */
    inline void         
    Reset               (const size_t* dim)                                      {

		for (size_t i = 0; i < INVALID_DIM; i++)
			_dim [i] = dim [i];

		Reset();

    }
    

    /**
     * @brief           Purge data and free RAM.
     */
    inline void         
    Clear               ()                                      {

    	for (size_t i = 0; i < INVALID_DIM; i++)
            _dim[i] = 1;

		_M.resize(0);

    }
    

    /**
     * @brief           Reset. i.e. Set all fields = T(0)
     */
    inline void         
    Zero               ()                                      {

		_M.assign (Size(), T(0)); 

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
    Matrix<T>&
    operator=           (const Matrix<T>& M);
    
    
    /**
     * @brief           Assignment operator. Sets all elements s.
     *
     * @param  s        The assigned scalar.
     */
    Matrix<T>           
    operator=           (const T& s);
    
    
    /**
     * @brief           Matrix product. i.e. this * M.
     *
     * @param  M        The factor.
     */
    Matrix<T>           
    operator->*         (const Matrix<T>& M) const;
   
    
    /**
     * @brief           Elementwise substruction of two matrices
     *
     * @param  M        Matrix substruent.
     */
    template <class S> Matrix<T>           
    operator-           (const Matrix<S>& M) const;
    
    
    /**
     * @brief           Elementwise substruction all elements by a scalar
     *
     * @param  s        Scalar substruent.
     */
    template <class S> Matrix<T>           
    operator-           (const S& s) const;
    
    
    /**
     * @brief           ELementwise substraction and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
	 * @return          Result
     */
    template <class S>  Matrix<T>           
    operator-=          (const Matrix<S>& M);
    
    
    /**
     * @brief           ELementwise substration with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator-=         (const S& s);
    
    
    /**
     * @brief           Unary minus (additive inverse)
	 *
	 * @return          Negation
     */
    Matrix<T>           
    operator-           () const;
    
    
    /**
     * @brief           Unary plus
	 *
	 * @return          Identity
     */
    Matrix<T>           
    operator+           () const;
    
    
    /**
     * @brief           Elementwise addition of two matrices
     *
     * @param  M        Matrix additive.
     */
    template <class S> Matrix<T>           
    operator+          (const Matrix<S>& M) const;
    
    
    /**
     * @brief           Elementwise addition iof all elements with a scalar
     *
     * @param  s        Scalar additive.
     */
    template <class S> Matrix<T>           
    operator+           (const S& s) const;
    
    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator+=          (const Matrix<S>& M);
    
    
    /**
     * @brief           ELementwise addition with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
	 * @return          Result
     */
    template <class S > Matrix<T>           
    operator+=          (const S& s);
    
    
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
	 * @param  M        The operand
	 * @return          Cross-section or zero
     */
    Matrix<T>           
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
    operator==          (const Matrix<T>& M) const;
    
    
    /**
     * @brief           Elementwise equality, result[i] = (this[i] != m[i]). i.e. this ~= m
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator!=          (const Matrix<T>& M) const;
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] >= m[i]). i.e. this >= m
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator>=          (const Matrix<T>& M) const;
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] <= m[i]). i.e. this <= m.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator<=          (const Matrix<T>& M) const;
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] > m[i]). i.e. this > m.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator>           (const Matrix<T>& M) const;
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] < m[i]). i.e. this < m.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>        
    operator<           (const Matrix<T>& M) const; 
    
    
    /**
     * @brief           Matrix comparison, result[i] = (m[i] || this[i] ? 1 : 0). i.e. this | m.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>           
    operator||          (const Matrix<T>& M) const;
    
    
    /**
     * @brief           Matrix comparison, result[i] = (m[i] && this[i] ? 1 : 0). i.e. this & m.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<bool>           
    operator&&          (const Matrix<T>& M) const;


    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Power.
	 * @return          Result
     */
    Matrix<T>           
    operator^           (const float& p) const;
    
    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Power.
	 * @return          Result
     */
    Matrix<T>           
    operator^=          (const float& p);
    

    /**
     * @brief           Elementwise multiplication. i.e. this .* M.
     *
     * @param  M        Factor matrix.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator*          (const Matrix<S> &M) const ;


    /**
     * @brief           Elementwise multiplication with a scalar. i.e. this * m.
	 *
	 * @param  s        Factor scalar
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator*          (const S& s) const ;

    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator*=         (const Matrix<S>& M);
    
    
    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this *= s.
     *
     * @param  s        Factor scalar.
	 * @return          Result
     */
    template <class S> Matrix<T>
    operator*=         (const S& s);
    
    
    /**
     * @brief           Elelemtwise division by M.
     *
     * @param  M        The divisor.
	 * @return          Result
     */
    template <class S>  Matrix<T>           
    operator/          (const Matrix<S>& M) const;

    
    /**
     * @brief           Elementwise division by scalar. i.e. this * m.
	 *
     * @param  s        The divisor.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator/           (const S& s) const;
    
    /**
     * @brief           ELementwise division and assignment operator. i.e. this = this ./ M.
     *
     * @param  M        Divisor matrix.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator/=         (const Matrix<S> &M);
    
    
    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Divisor scalar.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator/=         (const S& s);
    
    
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
	GetClassName        () const { return _name.c_str(); }


	/**
	 * @brief           Who are we?
	 *
	 * @return          Class name
	 */ 
	void 
	SetClassName        (const char* name) { _name = name; }


	/**
	 * @brief           Multiple subscripts from linear index (MATLAB like)
	 *
	 * @param  sub      Subscripts
	 * @return          Indices
	 */
	Matrix<size_t>
	Sub2Ind             (const Matrix<size_t>& sub) const;



	/**
	 * @brief           Multiple subscripts from linear index (MATLAB like)
	 *
	 * @param  ind      Indices
	 * @return          Subscripts
	 */
	Matrix<size_t>
	Ind2Sub             (const Matrix<size_t>& ind) const;



	/**
	 * @brief           Multiple subscripts from linear index (MATLAB like) (faster for 2D)
	 *
	 * @param  ind      Indices
	 * @return          Subscripts
	 */
	Matrix<size_t>
	Ind2Sub2D           (const Matrix<size_t>& ind) const;



	/**
	 * @brief           Multiple subscripts from linear index (MATLAB like) (faster for 3D)
	 *
	 * @param  ind      Indices
	 * @return          Subscripts
	 */
	Matrix<size_t>
	Ind2Sub3D           (const Matrix<size_t>& ind) const;



	/**
     * @brief           Matrix Product.
     *
	 * @param   M       The factor
	 * @param   transa  Transpose ('T') / Conjugate transpose ('C') the left matrix. Default: No transposition'N'
	 * @param   transb  Transpose ('T') / Conjugate transpose ('C') the right matrix. Default: No transposition 'N'
     * @return          Product of this and M.
     */
    Matrix<T>           
    prod                (const Matrix<T> &M, const char& transa = 'N', const char& transb = 'N') const;
    

    /**
     * @brief           Complex conjugate left and multiply with right.
     *
	 * @param   M       Factor
     * @return          Product of conj(this) and M.
     */
    Matrix<T>           
    prodt               (const Matrix<T> &M) const;
    

    /**
     * @brief           Complex conjugate right and multiply with right.
     *
	 * @param   M       Factor
     * @return          Product of conj(this) and M.
     */
    Matrix<T>           
    tprod               (const Matrix<T> &M) const;
    

	/**
     * @brief           Scalar product (complex: conjugate first vector) using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTC and XDOT
     *
     * @param  M        Factor
	 * @return          Scalar product
     */
	T
    dotc (const Matrix<T>& M) const;
    
	
	/**
     * @brief           Scalar product using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTU and XDOT
     *
     * @param  M        Factor
	 * @return          Scalar product
     */
	T
    dotu (const Matrix<T>& M) const;
    
	/**
     * @brief           Scalar product using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTU and XDOT
     *
     * @param  M        Factor
	 * @return          Scalar product
     */
	T
    dot (const Matrix<T>& M) const;
    
    /**
     * @brief           Transposition / Complex conjugation and transposition.
     *
     * @return          The transposed matrix
     */
    Matrix<T>           
    tr   ()             const;
    

	/**
	 * @brief           Number of elements
	 *
	 * @return          Size
	 */
	inline size_t
	Size() const {
		
		long size = 1;
		
		for (size_t i = 0; i < INVALID_DIM; i++)
			size *= _dim[i];
		
		return size;
		
	}


	/**
	 * @brief           Subscript of 1st dimension from index
	 *
	 * @param  ind      Index
	 * @return          Subscript of 1st dimension
	 */ 
	size_t 
	Ind2i               (const size_t& ind) const;


	/**
	 * @brief           Subscript of 2nd dimension from index
	 *
	 * @param  ind      Index
	 * @return          Subscript of 2nd dimension
	 */ 
	size_t 
	Ind2j               (const size_t& ind) const;


	/**
	 * @brief           Subscript of 3rd dimension from index
	 *
	 * @param  ind      Index
	 * @return          Subscript of 3rd dimension
	 */ 
	size_t 
	Ind2k               (const size_t& ind) const;


	/**
	 * @brief           Subscript of 4th dimension from index
	 *
	 * @param  ind      Index
	 * @return          Subscript of 4th dimension
	 */ 
	size_t 
	Ind2l               (const size_t& ind) const;


	/**
	 * @brief           Subscript of x-th dimension from index
	 *
	 * @param  ind      Index
	 * @param  dim      Dimension
	 * @return          Subscipt of x-th dimension
	 */ 
	size_t 
	Ind2x               (const size_t& ind, const size_t& dim) const;


    //@}


private:
    
    size_t              _dim[INVALID_DIM]; /// Dimensions
    float               _res[INVALID_DIM]; /// Resolutions
	std::vector<T>      _M;
	std::string         _name; 


	/**
	 * @brief           Adjust and resize for Syngo read
	 *
	 * @param  fname    Syngo MR meas file name
	 * @return          Success
	 */
	bool
	RSAdjust            (const std::string& fname);


    /**
     * @brief           Reset. i.e. reallocate and set all fields = T(0)
     */
    inline void         
    Reset               ()                                      {

		_M.resize(Size());
		Zero();

    }
    


    
};


template<class T> inline size_t
Matrix<T>::Ind2i  (const size_t& ind) const { 
	return (size_t) ind % _dim[0];                 
}


template<class T> inline size_t
Matrix<T>::Ind2j  (const size_t& ind) const { 
	return (size_t) floor (ind/_dim[0]) % (_dim[1]-1);
}


template<class T> inline size_t
Matrix<T>::Ind2k  (const size_t& ind) const { 
	return (size_t) floor (ind/(_dim[0]*_dim[1])) % (_dim[2]-1);
}


template<class T> inline size_t
Matrix<T>::Ind2l  (const size_t& ind) const { 
	return (size_t) floor (ind/(_dim[0]*_dim[1]*_dim[2])) % (_dim[3]-1);
}


template<class T> inline size_t
Matrix<T>::Ind2x (const size_t& ind, const size_t& dim) const { 
	
	size_t x = 1;

	for (size_t i = 1; i < dim+1; i++)
		x *= _dim[i-1]; 
	
	x = (size_t) floor((double)ind/(double)x) % (_dim[dim]);

	return (x);

}


template<class T> Matrix<size_t>
Matrix<T>::Ind2Sub2D (const Matrix<size_t>& inds) const {
	
	Matrix<T>      tmp = this->Squeeze();
	Matrix<size_t> subs (inds.Size(), 2);
	
	for(size_t i=0; i < subs.Width(); i++)
		for(size_t j=0; j < subs.Height() ; j++)
			subs(j,i) = Ind2x(inds(j), i);

	return subs; 

}


template <class T> inline Matrix<size_t>
Matrix<T>::Ind2Sub3D (const Matrix<size_t>& inds) const {
	
	Matrix <size_t> subs (inds.Size(), 3);
	
	for(size_t i=0; i < subs.Width(); i++)
		for(size_t j=0; j < subs.Height() ; j++)
			subs(j,i) = Ind2x(inds(j), i);

	return subs; 

}


template <class T> inline Matrix<size_t>
Matrix<T>::Sub2Ind  (const Matrix<size_t>& subs) const {

	size_t n = subs.Dim(0);

	Matrix<size_t> inds (n);

	/*for (int i = 0; i < n; i++)
	  inds[i] = */

	return subs; 
}


template <class T> inline Matrix<T> 
Matrix<T>::tr() const {

    Matrix<T> res (_dim);
	
	long tmp   = res.Dim(0);
	res.Dim(0) = res.Dim(1);
	res.Dim(1) = tmp;
	
    for (size_t i = 0; i < res.Dim(0); i++)
        for (size_t j = 0; j < res.Dim(1); j++)
			res.At(i,j) = cconj(At(j,i)); // Conjugate transpose
	
    return res;

}

#include "Lapack.hpp"

template <class T> Matrix<T> 
Matrix<T>::prodt (const Matrix<T> &M) const {
	
	return gemm (*this, M, 'C');
	
}


template <class T> Matrix<T> 
Matrix<T>::prod (const Matrix<T> &M, const char& transa, const char& transb) const {
	
	return gemm (*this, M, transa, transb);
	
}


template<class T> T 
Matrix<T>::dotc (const Matrix<T>& M) const {

	return DOTC (*this, M);
	
}

template<class T>  T 
Matrix<T>::dotu (const Matrix<T>& M) const {

	return DOTU (*this, M);
	
}

template<class T>  T 
Matrix<T>::dot (const Matrix<T>& M) const {
	
	return DOT  (*this, M);
	
}


template <class T> std::ostream& 
operator<< (std::ostream& os, Matrix<T>& M) {
	
	M.Print(os);
	return os;
	
}

template<> inline std::ostream&  
Matrix<size_t>::Print (std::ostream &os) const {
	
	for (size_t i = 0; i < _dim[COL]; i++) {
		for(size_t j = 0; j < _dim[LIN]; j++)
			printf ("%i ", (int)_M [i + j * _dim[COL]]);
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<short>::Print (std::ostream &os) const {
	
	for (size_t i = 0; i < _dim[COL]; i++) {
		for(size_t j = 0; j < _dim[LIN]; j++)
			printf ("%i ", _M [i + j * _dim[COL]]);
		printf("\n");
	}
	
	return os;

}


template<> inline std::ostream&  
Matrix<long>::Print (std::ostream &os) const {
	
	for (size_t i = 0; i < _dim[COL]; i++) {
		for(size_t j = 0; j < _dim[LIN]; j++)
			printf ("%li ", _M [i + j * _dim[COL]]);
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<double>::Print (std::ostream &os) const {
	
	for (size_t i = 0; i < _dim[COL]; i++) {
		for(size_t j = 0; j < _dim[LIN]; j++)
			printf ("%+.4f ", _M [i + j * _dim[COL]]);
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<float>::Print (std::ostream &os) const {
	
	for (size_t i = 0; i < _dim[COL]; i++) {
		for(size_t j = 0; j < _dim[LIN]; j++)
			printf ("%+.4f ", _M [i + j * _dim[COL]]);
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<cxfl>::Print (std::ostream& os) const {
	
	for (size_t i = 0; i < _dim[COL]; i++) {
		for(size_t j = 0; j < _dim[LIN]; j++)
			printf ("%+.4f+%+.4fi ", _M [i + j * _dim[COL]].real(), _M [i + j * _dim[COL]].imag());
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<cxdb>::Print (std::ostream& os) const {
	
	for (size_t i = 0; i < _dim[COL]; i++) {
		for(size_t j = 0; j < _dim[LIN]; j++)
			printf ("%+.4f+%+.4fi ", _M [i + j * _dim[COL]].real(), _M [i + j * _dim[COL]].imag());
		printf("\n");
	}
	
	return os;
	
}


template <class T> inline 
Matrix<T>::Matrix () {

    for (size_t i = 0; i < INVALID_DIM; i++) {
        _dim [i] = 1;
        _res [i] = 1.0;
	}

	_M.resize(Size());
		
	Zero();

	_name = "Matrix<T>";

}



template <class T> inline 
Matrix<T>::Matrix (const size_t& n) {

	_dim [COL] = n;
	_dim [LIN] = n;

    for (size_t i = 2; i < INVALID_DIM; i++)
        _dim [i] = 1;

	for (size_t i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;
	
	_M.resize(n*n, T(0));
	
}



template <class T> inline 
Matrix<T>::Matrix (const size_t& m, const size_t& n) {

	_dim [0] = m;
	_dim [1] = n;

    for (size_t i = 2; i < INVALID_DIM; i++) 
        _dim [i] = 1;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;

	_M.resize(m*n, T(0));

}



template <class T> inline 
Matrix<T>::Matrix (const size_t& m, const size_t& n, const size_t& k) {

	_dim [0] = m;
	_dim [1] = n;
	_dim [2] = k;
	
    for (size_t i = 3; i < INVALID_DIM; i++)
        _dim [i] = 1;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
        _res [i] = 1.0;
	
	_M.resize(m*n*k, T(0));
	
}



template <class T> inline 
Matrix<T>::Matrix (const size_t col, const size_t lin, const size_t cha, const size_t set, 
                   const size_t eco, const size_t phs, const size_t rep, const size_t seg, 
                   const size_t par, const size_t slc, const size_t ida, const size_t idb, 
                   const size_t idc, const size_t idd, const size_t ide, const size_t ave) {
	
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
	
	_M.resize(Size());
		
	Zero();

	
}



template <class T> inline 
Matrix<T>::Matrix (const size_t* dim) {
	
	for (size_t i = 0; i < INVALID_DIM; i++) {
		_dim[i] = dim[i];
        _res[i] = 1.0;
	}
	
	_M.resize(Size());
		
	Zero();

	
}


template <class T> inline 
Matrix<T>::Matrix (const size_t* d, const float* r) {
	
	for (size_t i = 0; i < INVALID_DIM; i++) {
		_dim[i] = d[i];
        _res[i] = r[i];
	}
	
	_M.resize(Size());
		
	Zero();

	
}


template <class T> inline 
Matrix<T>::Matrix (const Matrix<T> &M) {
	
	size_t i;

	for (i = 0; i < INVALID_DIM; i++) {
		_dim[i] = M.Dim(i);
		_res[i] = M.Res(i);
	}
	
	i = Size();
	_M.resize(i);
	
	memcpy (&_M[0], M.Data(), i * sizeof(T));
	
}



template <class T> inline 
Matrix<T>::~Matrix() {
    
#ifdef PARC_MODULE_NAME
    ICE_SET_FN ("Matrix<T>::~Matrix()")
    ICE_WARN   ("Freeing " << (float)Size() * sizeof(T) / 1024 << " kB of RAM.");
#endif

	_M.resize(0);
    
}



template <class T> inline Matrix<T>&
Matrix<T>::operator= (const Matrix<T>& M) {
    
	if (this != &M) {

		if (this->Size() != M.Size()) {

			try {
				_M.resize(M.Size());
			} catch (...) {
				printf ("Couldn't allocate %zu elements\n", M.Size()); 
			}
			
		}
		
		memcpy (_dim, M.Dim(), INVALID_DIM * sizeof(size_t));
		
		_M = M.Dat();
		
	}

    return *this;
	
}


template <class T> inline Matrix<T> 
Matrix<T>::operator= (const T& s) {
	
	_M.assign (Size(), s);
    return *this;
	
}


template <class T> inline Matrix<bool> 
Matrix<T>::operator== (const T& s) const {

    Matrix<bool> res(_dim);
	size_t i = Size();

	while (i--)
		res.Dat()[i] = (_M[i] == s);
	
    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator>= (const T& s) const {

    Matrix<bool> res(_dim);
	size_t i = Size();
    
	while (i--)
		res.Dat()[i] = (_M[i] >= s);

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator<= (const T& s) const {

    Matrix<bool> res(_dim);
	size_t i = Size();
		
	while (i--)
        res.Dat()[i] = (_M[i] <= s);

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator!= (const T& s) const {

    Matrix<bool> res(_dim);
	size_t i = Size();

	while (i--)
        res.Dat()[i] = (_M[i] != s);

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator< (const T& s) const {

    Matrix<bool> res(_dim);
	size_t i = Size();

	while (i--)
        res.Dat()[i] = (_M[i] < s);

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator> (const T& s) const {

    Matrix<bool> res(_dim);
	size_t i = Size();

	while (i--)
        res.Dat()[i] = (_M[i] > s);

    return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator->* (const Matrix<T> &M) const {

    return this->prod(M);

}


template <class T> inline Matrix<T> 
Matrix<T>::operator!() const {

    return this->tr();

}


template <class T> inline Matrix<T> 
Matrix<T>::operator& (const Matrix<bool>& M) const {

    for (size_t i = 0; i < INVALID_DIM; i++) 
        assert (_dim[i] == M.Dim()[i]);


    size_t k = 0, i = Size();

	while (i--)
        if (M[i])
            k++;
    
    Matrix<T> res(k, 1);
    
    k = 0; 
	i = Size();

	while (i--)
        if (M[i]) {
            res[k] = i;
            k++;
        }
	
    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator&& (const Matrix<T>& M) const {

	size_t i;

	for (i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<T> res(_dim);
	i = Size();

	while (i--)
        res[i] = (M[i] && _M[i]);

    return res;
	
}


template <class T> inline Matrix<bool> 
Matrix<T>::operator|| (const Matrix<T>& M) const {

	size_t i;

	for (i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);
	i = Size();

	while (i--)
        res.Dat()[i] = (_M[i] || M[i]);

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator== (const Matrix<T>& M) const {

	size_t i;

	for (i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);
	i = Size();

	while (i--)
        res.Dat()[i] = (_M[i] == M[i]);

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator>= (const Matrix<T>& M) const {

	size_t i;

	for (i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);
	i = Size();

	while (i--)
        res.Dat()[i] = (_M[i] >= M[i]) ? true : false;

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator<= (const Matrix<T>& M) const {

	size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);
	i = Size();

	while (i--)
        res.Dat()[i] = (_M[i] <= M[i]) ? true : false;

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator!= (const Matrix<T>& M) const {

	size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);
	i = Size();

	while (i--)
        res.Dat()[i] = (_M[i] != M[i]) ? true : false;

    return res;
	
}


template <class T> inline Matrix<bool> 
Matrix<T>::operator> (const Matrix<T>& M) const {
	
	size_t i;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));
	
    Matrix<bool> res(_dim);
	i = Size();

	while (i--)
        res.Dat()[i] = (_M[i] > M[i]) ? true : false;

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator< (const Matrix<T>& M) const {

	size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);
	i = Size();

	while (i--)
        res.Dat()[i] = (_M[i] < M[i]);
	
    return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator- () const {

	size_t i = Size();
    Matrix<T> res (_dim);

	while (i--)
        res[i] =- _M[i];
	
    return res;
	
}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator- (const Matrix<S> &M) const {
	
	size_t i;
	
	for (i=0; i < INVALID_DIM; i++)
		assert (Dim(i) == M.Dim(i));
	
    Matrix<T> res = *this;
	i = Size();

	while (i--)
		res[i] -= M[i];
	
	return res;
	
}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator- (const S& s) const {
	
    Matrix<T> res = *this;
	T t = T(s);
	size_t i = Size();
	
	while (i--)
		res[i] -= t;
	
	return res;
	
}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator+ (const Matrix<S> &M) const {

	size_t i;
	
	for (i=0; i < INVALID_DIM; i++)
		assert (Dim(i) == M.Dim(i));

    Matrix<T> res = *this;
	i = Size();

	while (i--)
		res[i] += M[i];
		
	return res;
	
}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator+ (const S& s) const {
	
    Matrix<T> res = *this;
	T t = T(s);
	size_t i = Size();
	
	while (i--)
		res[i] += t;
	
	return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator^ (const float& p) const {
    
	Matrix<T> res = *this;
	size_t i = Size();

	while (i--)
        res[i] = (p == 0) ? T(1) : pow(res[i],  p);

    return res;

}


template <class T> inline Matrix<T>
Matrix<T>::operator ^= (const float& p) {

	size_t i = Size();
		
	while (i--)
        _M[i] = (p == 0) ? T(1) : pow(_M[i],  p);
	
    return *this;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator* (const Matrix<S> &M) const {

	size_t i;

	for (i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<T> res = *this;
	i = Size();

	while(i--)
		res[i] *= M[i];

	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator* (const S& s) const {
    
    Matrix<T> res = *this;
	T t = T(s);
	size_t i = Size();

	while (i--)
		res[i] *= t;

	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator *= (const Matrix<S> &M) {
    
    size_t i;

	for (i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

	i = Size();
		
	while (i--)
		_M[i] *= M[i];
	
    return *this;

}


template <class T> template <class S> inline Matrix<T>
Matrix<T>::operator *= (const S& s) {
    
	T t = T(s);
	size_t i = Size();
	
	while (i--)
		_M[i] *= t;
	
    return *this;
	
}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator += (const Matrix<S> &M) {
	
    size_t i;
	
	for (i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

	i = Size();
	
	while(i--)
		_M[i] += M[i];
		
    return *this;

}


template <class T> template <class S> inline Matrix<T>
Matrix<T>::operator+= (const S& s) {
    
	T t = T(s);
	size_t i = Size();

	while (i--)
		_M[i] += t;
			
    return *this;
	
}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator-= (const Matrix<S>& M) {
	
    size_t i;
	
	for (i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));
	
	i = Size();
	
	while (i--)
		_M[i] -= M[i];
	
    return *this;
	
}


template <class T> template <class S> inline Matrix<T>
Matrix<T>::operator-= (const S& s) {
    
    size_t i = Size();
	T t = T(s);

	while (i--)
		_M[i] -= t;

    return *this;

}


template<class T> template<class S> inline Matrix<T> 
Matrix<T>::operator /= (const Matrix<S> &M) {
    
    size_t i;

	for (i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

	i = Size();

	while (i--)
		_M[i] /= M[i];

    return *this;

}


template <class T> template<class S> inline Matrix<T>
Matrix<T>::operator/= (const S& s) {
    
	T t = T(s);
	size_t i = Size();
		
	while (i--)
		_M[i] /= t;

    return *this;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator/ (const Matrix<S>& M) const {

	size_t i;

	for (i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

	i = Size();
    Matrix<T> res = *this;

	while(i--)
		res[i] = (M[i] != (T)0) ? _M[i] / M[i] : 0;

	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator/ (const S& s) const {
    
	assert (cabs(s) != 0.0);

	T t = T(s);
	size_t i = Size();
    Matrix<T> res = *this;
	
	while (i--)
		res[i] /= t;

	return res;

}


template <class T> inline T 
Matrix<T>::operator[]  (const size_t& p) const {
    
    assert(p <  Size());
    
    return _M[p];
    
}


template <class T> inline T&
Matrix<T>::operator[] (const size_t& p) {
    
    assert(p <  Size());
    
    return _M[p];
    
}


template<class T> template<class S> inline
Matrix<T>::operator Matrix<S> () const {

	Matrix<S> m (_dim);
	size_t i = Size();

	while (i--)
		m[i] = (S)_M[i];
		
	return m;

}

#ifdef PARC_MODULE_NAME

template <class T> long 
Matrix<T>::Import     (const IceAs* ias, const size_t pos) {
    
    ICE_SET_FN("Matrix<T>::Import(IceAs, long)")
        
    int  i    = 0;
    long size = 1;
    
    for (i = 0; i < INVALID_DIM; i++)
        size *= (ias->getLen(IceDim(i)) <= 1) ? 1 : ias->getLen(IceDim(i));
    
    T* data = (T*) ias->calcSplObjStartAddr() ;
    
    for (i = 0; i < size; i++, data++)
        _M[i+pos] = *data;
    
    return size;
    
}


template <class T> long 
Matrix<T>::Import(const IceAs* ias) {
    
    ICE_SET_FN("Matrix<T>::Import(IceAs)")
        
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


template <class T> long 
Matrix<T>::Export (IceAs* ias) const {
    
    ICE_SET_FN("Matrix<T>::Export(IceAs)")
		
    T* data = (T*) ias->calcSplObjStartAddr() ;
    
    for (int i = 0; i < Size(); i++, data++)
        *data = _M[i];
    
    return Size();
    
}


template <class T> long
Matrix<T>::Export (IceAs* ias, const size_t pos) const {

    ICE_SET_FN("Matrix<T>::Export(IceAs, long)")
        
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

#endif // __MATRIX_H__
