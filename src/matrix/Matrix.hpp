/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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
    # define PI  3.1415926535897931159979634685441851615906
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


enum io_strategy {

	HDF5,
	MATLAB,
	NIFTI,
	SYNGO,
	CDF,
	PRIMITIVE

};


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


	/**
	 * @brief           Create nxn identity matrix 
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn identity
	 */
	static Matrix<T> 
	Id                  (const size_t n);
	
	
	/**
	 * @brief           Create nxn matrix initialised with T(1.0)
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn ones
	 */
	static Matrix<T> 
	Ones                (const size_t n);
	
	
	/**
	 * @brief           Create mxn matrix initialised with T(1.0)
	 *
	 * @param  m        Rows
	 * @param  n        Cols
	 * @return          mxn ones
	 */
	static Matrix<T> 
	Ones                (const size_t m, const size_t n);
	
	
	/**
	 * @brief           Create mxnxl matrix initialised with T(1.0)
	 *
	 * @param  m        Rows
	 * @param  n        Cols
	 * @param  l        Slices
	 * @return          nxn ones
	 */
	static Matrix<T> 
	Ones                (const size_t m, const size_t n, const size_t l);
	
	
	/**
	 * @brief           Create nxn matrix initialised with T(0.0)
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn zeros
	 */
	static Matrix<T> 
	Zeros               (const size_t n);
	
	
	/**
	 * @brief           Create mxn matrix initialised with T(0.0)
	 *
	 * @param  m        Rows
	 * @param  n        Cols
	 * @return          nxn zeros
	 */
	static Matrix<T> 
	Zeros               (const size_t m, const size_t n);
	
	
	/**
	 * @brief           Create nxn matrix initialised with T(0.0)
	 *
	 * @param  m        Rows
	 * @param  n        Cols
	 * @param  l        Slices
	 * @return          nxn zeros
	 */
	static Matrix<T> 
	Zeros               (const size_t m, const size_t n, const size_t l);

	
	/**
	 * @brief           Create 2D nxn Shepp-Logan phantom.<br/>
	 *                  Shepp et al. The Fourier reconstruction of a head section. IEEE TNS. 1974; 21: 21-43
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn zeros
	 */
	static Matrix<T> 
	Phantom2D           (const size_t n);
	
	
	/**
	 * @brief           Create 3D nxn Shepp-Logan phantom.<br/>
	 *                  Koay et al. Three dimensional analytical magnetic resonance imaging phantom in the Fourier domain. MRM. 2007; 58: 430-436
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn zeros
	 */
	static Matrix<T> 
	Phantom3D           (const size_t n);
	
	
	/**
	 * @brief           Create 2D nxn random matrix.<br/>
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn zeros
	 */
	static Matrix<T> 
	Random2D            (const size_t n);
	
	
	/**
	 * @brief           Create 2D nxn random matrix.<br/>
	 *
	 * @param  n        Side length of matrix
	 * @return          nxn zeros
	 */
	static Matrix<T> 
	Random3D            (const size_t n);
	
	
	/**
	 * @brief           Create circle of ones in 2D square Matrix
	 *
	 * @param   p       Parameter array. Must hold r, x0, y0, intensity.
	 * @param   n       Size of square
	 * @return          Matrix including ellipsoid
	 */
	static Matrix<T> 
	Circle              (const float* p, const size_t n);
	
	
	/**
	 * @brief           Create sphere of ones in 3D cubic Matrix
	 *
	 * @param   p       Parameter array. Must hold r, x0, y0, z0, size_tensity.
	 * @param   n       Size of cube
	 * @return          Matrix including ellipsoid
	 */
	static Matrix<T> 
	Sphere              (const float* p, const size_t n);
	
	
	/**
	 * @brief           Create ellipsoid of ones in 3D cubic matrix
	 *
	 * @param   p       Parameter array. Must hold a, b, c, x0, y0, z0, phi, psi, theta, intensity.
	 * @param   n       Size of cube
	 * @param   v       Value of voxels inside, default: 1
	 * @return          Matrix including ellipsoid
	 */
	static Matrix<T> 
	Ellipsoid           (const float* p, const size_t n, const T v = T(1.0));
	
	
	/**
	 * @brief           Create ellipse of ones in 2D square matrix
	 *
	 * @param  p        Parameter array. Must hold a, b, x0, y0, phi, intensity.
	 * @param  n        Size of square
	 * @param  v        Value of voxels inside, default: 1
	 * @return          Matrix including ellipsoid
	 */
	static Matrix<T> 
	Ellipse             (const float* p, const size_t n, const T v = T(1.0));
	
	
	/**
	 * @brief           MATLAB-like linspace
	 *
	 * @param  start    Start of range
	 * @param  end      End of range
	 * @param  n        Number of samples
	 * @return          Vector of values
	 */
	static Matrix<T> 
	LinSpace           (const T& start, const T& end, const size_t& n);
	
	
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
    operator()          (const size_t& p) const;

    
    /**
     * @brief           Get value of pth element of repository.
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    T&                 
    operator()          (const size_t& p) ;

    
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
        return _M[col + _dim[COL]*lin ];
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
    operator()           (const size_t& col, const size_t& lin, const size_t& slc) const {
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
    inline T&                 
    operator()           (const size_t& col, const size_t& lin, const size_t& slc) {
           return _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc];
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

		for (size_t i = 0; i < res.Size(); i++)
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
    operator->*         (Matrix<T>& M);
   
    
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
    Matrix<T>           
    operator||          (const Matrix<T>& M) const;
    
    
    /**
     * @brief           Matrix comparison, result[i] = (m[i] && this[i] ? 1 : 0). i.e. this & m.
     *
     * @param  M        Comparing matrix.
	 * @return          Hit list
     */
    Matrix<T>           
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
     * @name Data manipulation functions.
     *       Data manipulation functions.
     */
    
    //@{
    
    /**
     * @brief           Fill with random data
     */
    inline void         
	Random             ();    


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
     * @brief           Greatest element of the matrix. i.e. max(M0).
     *
     * @return          The scalar value of the greatest element.
     */
    T                   
    Max                 () const ;    
    
    
    /**
     * @brief           Smallest element of the matrix. i.e. min(M0).
     *
     * @return          The scalar value of the smallest element.
     */
    T                   
    Min                 () const ;
    
    
	/**
	 * @brief           MeshGrid
	 *
	 * @param  dims     Grid indices (i.e. m(0,0)=minx m(0,1)=maxx, etc)
	 * @return          Mesh grid.<br/> [X,Y,Z] = meshgrid (.,.,.) corresponds to:<br/> 
	 *                  Matrix<size_t> dims (3,2); <br/>
	 *                  dims (0,0)=0; dims(0,1)=2; dims(1,0)=0; dims(1,1)=3; dims(2,0)=0; dims(2,1)=4;
	 *                  Matrix<size_t> mg = Matrix<T>::MeshGrid(dims);<br/>
	 *                  Delivers mg O(3,4,5,3). Highest dimension is [X,Y,Z];
	 */
	static Matrix<size_t>
	MeshGrid            (const Matrix<size_t>& dims);


    /**
     * @brief           Get maximum absolute value in the matrix
     *
     * @return          Maximum value
     */
    T                   
    Maxabs              () const ;
    
    
    /**
     * @brief           Get minimum absolute value in the matrix
     *
     * @return          Maximum value
     */
    T                   
    Minabs              () const ;
    

	/**
     * @brief           Matrix Product.
     *
	 * @param   M       The factor
	 * @param   transa  Transpose ('T') / Conjugate transpose ('C') the left matrix. Default: No transposition'N'
	 * @param   transb  Transpose ('T') / Conjugate transpose ('C') the right matrix. Default: No transposition 'N'
     * @return          Product of this and M.
     */
    Matrix<T>           
    prod                (Matrix<T> &M, const char transa = 'N', const char transb = 'N');
    
    /**
     * @brief           Complex conjugate left and multiply with right.
     *
	 * @param   M       Factor
     * @return          Product of conj(this) and M.
     */
    Matrix<T>           
    prodt               (Matrix<T> &M);
    

	/**
     * @brief           Scalar product (complex: conjugate first vector) using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTC and XDOT
     *
     * @param  M        Factor
	 * @return          Scalar product
     */
	T
    dotc (Matrix<T>& M);
    
	
	/**
     * @brief           Scalar product using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTU and XDOT
     *
     * @param  M        Factor
	 * @return          Scalar product
     */
	T
    dotu (Matrix<T>& M);
    
	/**
     * @brief           Scalar product using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTU and XDOT
     *
     * @param  M        Factor
	 * @return          Scalar product
     */
	T
    dot (Matrix<T>& M);
    
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


template <> inline short 
Matrix<short>::Max() const {
	
    short max = _M[0];
	
    for (size_t i = 0; i < Size(); i++)
        if (_M[i] > max)
            max = _M[i];
	
    return max;
	
}


template <> inline double 
Matrix<double>::Max() const {
	
    short max = _M[0];
	
    for (size_t i = 1; i < Size(); i++)
        if (_M[i] > max)
            max = _M[i];
	
    return max;
	
}


template <> inline cxfl
Matrix<cxfl>::Max() const {
		
	cxfl   max = cxfl(0.0,0.0);
	float tmp =  0.0;
	
	for (size_t i = 0; i < Size(); i++) {
		float abs = sqrt(_M[0].real()*_M[0].real() + _M[0].imag()*_M[0].imag());
		if (abs > tmp) {
			tmp = abs;
			max = _M[i];
		}
	}
	return max;
	
}


template <class T> inline T
Matrix<T>::Maxabs() const {

    T max = abs(_M[0]);

    for (size_t i = 0; i < Size(); i++)
        if (abs(_M[i]) > abs(max))
            max = abs(_M[i]);

    return cabs(max);

}


template <class T> inline T
Matrix<T>::Min() const {

    T min = _M[0];

    for (size_t i = 0; i < Size(); i++)
        if (_M[i] < min)
            min = _M[i];

    return min;

}


template <class T> inline T  
Matrix<T>::Minabs() const {

    T old = fabs(_M[0]);

    for (size_t i = 0; i < Size(); i++)
        if (fabs(_M[i]) < old)
            old = cabs(_M[i]);

    return old;

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

template<> inline void 
Matrix<cxfl>::Random () {
	
	srand (time(NULL));

	for (size_t i = 0; i < Size(); i++)
		_M[i] = cxfl ((float) rand() / (float) RAND_MAX*2-1, (float) rand() / (float) RAND_MAX*2-1);
	
}
    
template<> inline void 
Matrix<cxdb>::Random () {
	
	srand (time(NULL));

	for (size_t i = 0; i < Size(); i++)
		_M[i] = cxdb ((float) rand() / (float) RAND_MAX*2-1, (float) rand() / (float) RAND_MAX*2-1);
	
}
    
template<> inline void 
Matrix<double>::Random () {

	srand (time(NULL));

	for (size_t i = 0; i < Size(); i++)
		_M[i] = (double) rand() / (double) RAND_MAX*2-1;

}
    
template<> inline void 
Matrix<short>::Random () {

	srand (time(NULL));

	for (size_t i = 0; i < Size(); i++)
		_M[i] = (short) 12 * (double)rand() / (double)RAND_MAX*2-1;

}


#include "Lapack.hpp"

template <class T> Matrix<T> 
Matrix<T>::prodt (Matrix<T> &M) {
	
	return GEMM (*this, M, 'C');
	
}


template <class T> Matrix<T> 
Matrix<T>::prod (Matrix<T> &M, const char transa, const char transb) {
	
	return GEMM (*this, M, transa, transb);
	
}


template<class T>  T 
Matrix<T>::dotc (Matrix<T>& M)  {

	return DOTC (*this, M);
	
}

template<class T>  T 
Matrix<T>::dotu (Matrix<T>& M)  {

	return DOTU (*this, M);
	
}

template<class T>  T 
Matrix<T>::dot (Matrix<T>& M)  {
	
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


template <class T> Matrix<T>
Matrix<T>::LinSpace (const T& start, const T& end, const size_t& n) {
	
	assert (n > 1);
	
	Matrix<T> res (n, 1);
	T gap;

	gap      = T(end-start) / T(n-1);
	
	res[0]   = start;
	res[n-1] = end;
	
	for (int i = 1; i < n-1; i++)
		res[i] = res[i-1] + gap;
	
	return res;
	
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
	
	for (size_t i = 0; i < INVALID_DIM; i++) {
		_dim[i] = M.Dim(i);
		_res[i] = M.Res(i);
	}
	   
	_M.resize(Size());
	
#pragma omp parallel default (shared)
	{
#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			_M[i] = M[i];
	}
			//	memcpy (&_M[0], M.Data(), Size() * sizeof(T));
	
}



template <class T> inline 
Matrix<T>::~Matrix() {
    
#ifdef PARC_MODULE_NAME
    ICE_SET_FN ("Matrix<T>::~Matrix()")
    ICE_WARN   ("Freeing " << (float)Size() * sizeof(T) / 1024 << " kB of RAM.");
#endif

	_M.resize(0);
    
}



template <class T> inline 
Matrix<T> Matrix<T>::Id (const size_t n) {

 	Matrix<T> M (n);

 	for (size_t i = 0; i < n; i++)
 		M[i*n+i] = T(1.0);

 	return M;

}



template <class T> inline 
Matrix<T> Matrix<T>::Ones (const size_t m, const size_t n, const size_t l) {

 	Matrix<T> M (m,n,l);

 	for (size_t i = 0; i < M.Size(); i++)
 		M[i] = T(1.0);

 	return M;

}



template <class T> inline 
Matrix<T> Matrix<T>::Ones (const size_t m, const size_t n) {

 	Matrix<T> M (m,n);

 	for (size_t i = 0; i < M.Size(); i++)
 		M[i] = T(1.0);

 	return M;

}



template <class T> inline 
Matrix<T> Matrix<T>::Ones (const size_t n) {

 	return Ones(n,n);

}



template <class T> inline 
Matrix<T> Matrix<T>::Zeros (const size_t n, const size_t m, const size_t l) {

 	Matrix<T> M (m,n,l);

 	for (size_t i = 0; i < M.Size(); i++)
 		M[i] = T(0.0);

 	return M;

}



template <class T> inline 
Matrix<T> Matrix<T>::Zeros (const size_t n, const size_t m) {

 	Matrix<T> M (m,n);

 	for (size_t i = 0; i < M.Size(); i++)
 		M[i] = T(0.0);

 	return M;

}



template <class T> inline 
Matrix<T> Matrix<T>::Zeros (const size_t n) {

 	return Zeros(n,n);

}



template <class T> inline 
Matrix<T> Matrix<T>::Circle (const float* p, const size_t n) {

	Matrix<T> res = Matrix<T>::Zeros(n);

	float m[2];
	float rad;

	rad = p[0] * float(n) / 2.0;

	m[0] = (1.0 - p[1]) * float(n) / 2.0;
	m[1] = (1.0 - p[2]) * float(n) / 2.0;

	for (size_t r = 0; r < res.Dim(1); r++)
		for (size_t c = 0; c < res.Dim(0); c++)
			res(c,r) = ( pow(((float)c-m[0])/rad, 2.0 ) + pow(((float)r-m[0])/rad, 2.0) <= 1.0) ? T(1.0) : T(0.0);

	return res;

}



template <class T> inline 
Matrix<T> Matrix<T>::Sphere (const float* p, const size_t n) {

	Matrix<T> res = Matrix<T>::Zeros(n,n,n);

	float m[3];
	float rad;

	rad = p[0] * float(n) / 2.0;

	m[0] = (1.0 - p[1]) * float(n) / 2.0;
	m[1] = (1.0 - p[2]) * float(n) / 2.0;
	m[2] = (1.0 - p[3]) * float(n) / 2.0;

	for (size_t s = 0; s < res.Dim(2); s++)
		for (size_t r = 0; r < res.Dim(1); r++)
			for (size_t c = 0; c < res.Dim(0); c++)
				res(c,r) = ( pow (((float)c-m[0])/rad, 2.0) + pow (((float)r-m[1])/rad, 2.0) + pow (((float)s-m[2])/rad, 2.0) <= 1.0) ? T(1.0) : T(0.0);

	return res;

}



template <class T> inline 
Matrix<T> Matrix<T>::Ellipse (const float* p, const size_t n, const T v) {

	Matrix<T> res = Matrix<T>::Zeros(n);

	float m[2];
	float a[2];

	a[0] = p[0] * float(n) / 2.0;
	a[1] = p[1] * float(n) / 2.0;

	m[0] = (1.0 - p[2]) * float(n) / 2.0;
	m[1] = (1.0 - p[3]) * float(n) / 2.0;

	float cosp = cos(p[4]);
	float sinp = sin(p[4]);
	
#pragma omp parallel default (shared) 
	{
		
		size_t tid      = omp_get_thread_num();
		size_t chunk    = n / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk) 
		
	for (size_t r = 0; r < n; r++)
		for (size_t c = 0; c < n; c++)
			res(c,r) = (pow( (((float)c-m[1])*cosp+((float)r-m[0])*sinp)/a[1], 2.0 ) + 
						pow( (((float)r-m[0])*cosp-((float)c-m[1])*sinp)/a[0], 2.0) <= 1.0) ? v : T(0.0);

	}

	return res;

}



template <class T> inline 
Matrix<T> Matrix<T>::Ellipsoid (const float* p, const size_t n, const T v) {

	Matrix<T> res = Matrix<T>::Zeros(n,n,n);

	float m[3];
	float a[3];
	float d;

	a[0] = p[0] * float(n) / 2.0;
	a[1] = p[1] * float(n) / 2.0;
	a[2] = p[2] * float(n) / 2.0;

	m[0] = (1.0 - p[3]) * float(n) / 2.0;
	m[1] = (1.0 - p[4]) * float(n) / 2.0;
	m[2] = (1.0 - p[5]) * float(n) / 2.0;

	float cosp = cos(p[6]);
	float sinp = sin(p[6]);
	
#pragma omp parallel default (shared) 
	{
		
		size_t tid      = omp_get_thread_num();
		size_t chunk    = n / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk) 
		
		for (size_t s = 0; s < n; s++)
			for (size_t r = 0; r < n; r++)
				for (size_t c = 0; c < n; c++)
					res(c,r,s) = ( pow( (((float)c-m[1])*cosp+((float)r-m[0])*sinp)/a[1], 2.0) + 
								   pow( (((float)r-m[0])*cosp-((float)c-m[1])*sinp)/a[0], 2.0) +
								   pow( ((float)s-m[2])/a[2], 2.0) <= 1.0) ? v : T(0.0);
		
	}

	return res;

}



template <class T> inline 
Matrix<T> Matrix<T>::Phantom2D (const size_t n) {

	const size_t ne = 10; // Number of ellipses
	const size_t np = 5;  // Number of geometrical parameters

	float p[ne][np] = {
		{ 0.6900, 0.9200,  0.00,  0.0000,  0.0 },
		{ 0.6624, 0.8740,  0.00, -0.0184,  0.0 },
        { 0.1100, 0.3100, -0.22,  0.0000, -0.3 },
		{ 0.1600, 0.4100,  0.22,  0.0000,  0.3 },
		{ 0.2100, 0.2500,  0.00,  0.3500,  0.0 },
		{ 0.0460, 0.0460,  0.00,  0.1000,  0.0 },
		{ 0.0460, 0.0460,  0.00, -0.1000,  0.0 },
		{ 0.0460, 0.0230,  0.08, -0.6050,  0.0 },
		{ 0.0230, 0.0230,  0.00, -0.6060,  0.0 },
		{ 0.0230, 0.0460, -0.06, -0.6050,  0.0 }
	};

	// Size_Tensities
	T v[ne] = {T(1.0), T(-0.8), T(-0.2), T(-0.2), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1)};

	// Empty matrix
	Matrix<T> res = Matrix<T>::Zeros(n);
	Matrix<T>        e;

	for (size_t i = 0; i < ne; i++) {
		e    = Matrix<T>::Ellipse (p[i], n, v[i]);
		res += e;
	}

	return res;

}


template <class T> inline 
Matrix<T> Matrix<T>::Random2D (const size_t n) {

	Matrix<T> res (n);
	res.Random();

	return res;

}



template <class T> inline 
Matrix<T> Matrix<T>::Random3D (const size_t n) {

	Matrix<T> res (n, n, n);
	res.Random();

	return res;

}



template <class T> inline 
Matrix<T> Matrix<T>::Phantom3D (const size_t n) {

	const size_t ne = 10; // Number of ellipses
	const size_t np =  9; // Number of geometrical parameters

	float p[ne][np] = {
		{ 0.690, 0.920, 0.900,  0.00,  0.000,  0.000,  0.0, 0.0, 0.0 },
        { 0.662, 0.874, 0.880,  0.00,  0.000,  0.000,  0.0, 0.0, 0.0 },
        { 0.110, 0.310, 0.220, -0.22,  0.000, -0.250, -0.3, 0.0, 0.0 },
        { 0.160, 0.410, 0.210,  0.22,  0.000, -0.250,  0.3, 0.0, 0.0 },
        { 0.210, 0.250, 0.500,  0.00,  0.350, -0.250,  0.0, 0.0, 0.0 },
        { 0.046, 0.046, 0.046,  0.00,  0.100, -0.250,  0.0, 0.0, 0.0 },
        { 0.046, 0.023, 0.020,  0.08, -0.650, -0.250,  0.0, 0.0, 0.0 },
        { 0.046, 0.023, 0.020,  0.06, -0.650, -0.250,  0.0, 0.0, 0.0 },
        { 0.056, 0.040, 0.100, -0.06, -0.105,  0.625,  0.0, 0.0, 0.0 },
        { 0.056, 0.056, 0.100,  0.00,  0.100,  0.625,  0.0, 0.0, 0.0 }
	};

	T v[ne] = {2.0, -0.8, -0.2, -0.2, 0.2, 0.2, 0.1, 0.1, 0.2, -0.2};

	Matrix<T> res = Matrix<T>::Zeros(n,n,n);
	Matrix<T> e;
	
	for (size_t i = 0; i < ne; i++) {
		e    = Matrix<T>::Ellipsoid (p[i], n, v[i]);
		res += e;
	}

	return res;

}



template<class T> Matrix<size_t>
Matrix<T>::MeshGrid (const Matrix<size_t>& d) {

	size_t side [3];

	side[0] = d(0,1) - d(0,0) + 1;
	side[1] = d(1,1) - d(1,0) + 1;
	side[2] = d(2,1) - d(2,0) + 1;

    Matrix<size_t> mg (side[1], side[0], side[2], 3);

	for (size_t s = 0; s < side[2]; s++)
		for (size_t l = 0; l < side[0]; l++)
			for (size_t c = 0; c < side[1]; c++) {
				mg(c,l,s,0) = l + d(0,0);
				mg(c,l,s,1) = c + d(0,1);
				mg(c,l,s,2) = s + d(0,2);
			}
	
	return mg;

}


template <class T> inline Matrix<T>&
Matrix<T>::operator= (const Matrix<T>& M) {
    
	if (this->Size() != M.Size())
		_M.resize(M.Size());
	
	memcpy (_dim, M.Dim(), INVALID_DIM * sizeof(size_t));

	_M = M.Dat();

    return *this;
	
}


template <class T> inline Matrix<T> 
Matrix<T>::operator= (const T& s) {
	
	_M.assign (Size(), s);
    return *this;
	
}


template <class T> inline Matrix<bool> 
Matrix<T>::operator== (const T& s) const   {

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			res[i] = (_M[i] == s);
		
	}
	
    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator>= (const T& s) const {

    Matrix<bool> res(_dim);
    
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			res[i] = (_M[i] >= s);

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator<= (const T& s) const {

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] <= s);
    

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator!= (const T& s) const {

    Matrix<bool> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] != s);

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator< (const T& s) const {

    Matrix<bool> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] < s);

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator> (const T& s) const {

    Matrix<bool> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] > s);

    return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator->* (Matrix<T> &M) {

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


    size_t k = 0, i;
    for (i = 0; i < Size(); i++)
        if (M[i])
            k++;
    
    Matrix<T> res(k, 1);
    
    k = 0;
    for (i = 0; i < Size(); i++)
        if (M[i]) {
            res[k] = i;
            k++;
        }

    return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator&& (const Matrix<T>& M) const {

    assert(M.Size() == Size());

    Matrix<T> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (M[i] && _M[i]);

    return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator|| (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] || M[i]) ? true : false;

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator== (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] == M[i]) ? true : false;
	
	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator>= (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] >= M[i]) ? true : false;

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator<= (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] <= M[i]) ? true : false;

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator!= (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] != M[i]) ? true : false;

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator> (const Matrix<T>& M) const {
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] > M[i]) ? true : false;

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator< (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] < M[i]);

	}

    return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator- () const {

    Matrix<T> res (this->Dim());

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] =- _M[i];

	}

    return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator- (const Matrix<S> &M) const {

	for (size_t i=0; i < INVALID_DIM; i++)
		assert (Dim(i) == M.Dim(i));
	
    Matrix<T> res = *this;
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			res[i] -= M[i];
		
	}
	
	return res;
	
}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator- (const S& s) const {
	
    Matrix<T> res = *this;
	T t = T(s);
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			res[i] -= t;
		
	}
	
	return res;
	
}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator+ (const Matrix<S> &M) const {

	for (size_t i=0; i < INVALID_DIM; i++)
		assert (Dim(i) == M.Dim(i));

    Matrix<T> res = *this;

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			res[i] += M[i];
		
	}
	
	return res;
	
}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator+ (const S& s) const {
	
    Matrix<T> res = *this;
	T t = T(s);
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			res[i] += t;
		
	}
	
	return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator^ (const float& p) const {
    
	Matrix<T> res = *this;
    
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
	for (size_t i = 0; i < Size(); i++)
        res [i] = (p == 0) ? res[i] = 1 : pow(res[i],  p);

	}

    return res;

}


template <class T> inline Matrix<T>
Matrix<T>::operator ^= (const float& p) {
    
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
	for (size_t i = 0; i < Size(); i++)
        _M[i] = (p == 0) ? T(1) : pow(_M[i],  p);
	
	}

    return *this;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator* (const Matrix<S> &M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<T> res = *this;

#pragma omp parallel default (shared) 
	{

#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			res[i] *= M[i];

	}
	
	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator* (const S& s) const {
    
    Matrix<T> res = *this;
	T t = T(s);

#pragma omp parallel default (shared) 
	{

#pragma omp for
	for (size_t i = 0; i < Size(); i++)
		res[i] *= t;
	
	}

	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator *= (const Matrix<S> &M) {
    
    size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] *= M[i];
		
	}

    return *this;

}


template <class T> template <class S> inline Matrix<T>
Matrix<T>::operator *= (const S& s) {
    
	T t = T(s);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] *= t;
		
	}

    return *this;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator += (const Matrix<S> &M) {
    
    size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] += M[i];
		
	}

    return *this;

}


template <class T> template <class S> inline Matrix<T>
Matrix<T>::operator+= (const S& s) {
    
	T t = T(s);

#pragma omp parallel default (shared) 
		{
			
#pragma omp for
			
			for (size_t i = 0; i < Size(); i++)
				_M[i] += t;
			
		}
		
    return *this;
	
}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator-= (const Matrix<S>& M) {
	
    size_t i;
	
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


template <class T> template <class S> inline Matrix<T>
Matrix<T>::operator-= (const S& s) {
    
    size_t i;
	T t = T(s);

#pragma omp parallel default (shared) 
		{
			
#pragma omp for
			
			for (size_t i = 0; i < Size(); i++)
				_M[i] -= t;
			
		}

    return *this;

}


template<class T> template<class S> inline Matrix<T> 
Matrix<T>::operator /= (const Matrix<S> &M) {
    
    size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] /= M[i];
		
	}

    return *this;

}


template <class T> template<class S> inline Matrix<T>
Matrix<T>::operator/= (const S& s) {
    
	T t = T(s);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] /= s;
		
	}

    return *this;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator/ (const Matrix<S>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<T> res = *this;

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			(M[i] != (T)0) ? res[i] = _M[i] / M[i] : 0;

	}

	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator/ (const S& s) const {
    
	assert (cabs(s) != 0.0);

	T t = T(s);

    Matrix<T> res = *this;

#pragma omp parallel default (shared) 
	{
		
#pragma omp for 
		
	for (size_t i = 0; i < Size(); i++)
		res[i] /= t;

	}

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


template <class T> inline T 
Matrix<T>::operator() (const size_t& a) const {

    assert(a <  Size());

    return _M[a];

}


template<class T> template<class S> inline
Matrix<T>::operator Matrix<S> () const {

	Matrix<S> m (_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (int i = 0; i < this->Size(); i++)
			m[i] = (S)_M[i];
		
	}
	
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
