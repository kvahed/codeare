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

/**
 * @brief Interpolation methods
 */
enum InterpMethod {
	LINEAR, BSPLINE
};

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
 * Return absolute value.
 */
# define ABS(A) (A > 0 ? A : -A)


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
    #define GAMMA 4.2576e7
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
	 * @brief           Construct square 2D matrix
	 *
	 * @param  n        Rows & Columns
	 */
    inline              
    Matrix              (const size_t n) ;
    
    
    /**
	 * @brief           Construct 2D matrix
	 *
	 * @param  m        Rows
	 * @param  n        Columns
	 */
	inline 
	Matrix              (const size_t m, const size_t n);
	
    
    /**
	 * @brief           Construct 3D volume
	 *
	 * @param  m        Rows
	 * @param  n        Columns
	 * @param  k        Slices
	 */
	inline 
	Matrix              (const size_t m, const size_t n, const size_t k);
	

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
	 * @brief           Create ellipse of ones in 2D square matrix
	 *
	 * @param  p        Parameter array. Must hold a, b, x0, y0, phi, intensity.
	 * @param  n        Size of square
	 * @param  v        Value of voxels inside, default: 1
	 * @return          Matrix including ellipsoid
	 */
	static Matrix<T> 
	LinSpace           (const T& start, const T& space, const T& end);
	
	
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
    operator[]          (const size_t p)                             const;
    
    
    /**
     * @brief           Reference to pth element from data repository.
     *
     * @param  p        Requested position.
     * @return          Reference to _M[p].
     */
    T                   
    &operator[]         (const size_t p)                              ;

    
    /**
     * @brief           Get pointer to data
     *  
     * @return          Data 
     */
    inline T*            
    Data                ()  {
        return &(_M.at(0));
	}

    
    /**
     * @brief           Get element at position 
     *  
     * @param  pos      Position
     * @return          Value at _M[pos]
     */
    inline T            
    At                  (const size_t pos)  const {

        return _M[pos];

    }

    

    /**
     * @brief            Reference to value at position
     *  
     * @param  pos       Position
     * @return           Reference to _M[pos]
     */
    inline T&           
    At                  (const size_t pos) {

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
    At                  (const size_t col, const size_t lin) const {

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
    At                  (size_t col, size_t lin) {

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
    At                   (size_t col, size_t lin, size_t slc)  const {

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
    At                   (size_t col, size_t lin, size_t slc) {

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
    At                   (const size_t col, 
						  const size_t lin, 
						  const size_t cha,
						  const size_t set,
						  const size_t eco,
						  const size_t phs = 0,
						  const size_t rep = 0,
						  const size_t seg = 0,
						  const size_t par = 0,
						  const size_t slc = 0,
						  const size_t ida = 0,
						  const size_t idb = 0,
						  const size_t idc = 0,
						  const size_t idd = 0,
						  const size_t ide = 0,
						  const size_t ave = 0) const {
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
    At                   (const size_t col, 
						  const size_t lin, 
						  const size_t cha,
						  const size_t set,
						  const size_t eco = 0,
						  const size_t phs = 0,
						  const size_t rep = 0,
						  const size_t seg = 0,
						  const size_t par = 0,
						  const size_t slc = 0,
						  const size_t ida = 0,
						  const size_t idb = 0,
						  const size_t idc = 0,
						  const size_t idd = 0,
						  const size_t ide = 0,
						  const size_t ave = 0) {

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
	Reshape             (const size_t col, 
						 const size_t lin, 
						 const size_t cha = 1,
						 const size_t set = 1,
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
						 const size_t ave = 1) const {
		
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
	Reshape             (const size_t col, 
						 const size_t lin, 
						 const size_t cha = 1,
						 const size_t set = 1,
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
						 const size_t ave = 1) {
		
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
    operator()          (const size_t p) const;

    
    /**
     * @brief           Get value of pth element of repository.
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    T&                 
    operator()          (const size_t p) ;

    
    /**
	 * @brief           Get value in slice
	 *
	 * @param  col      Column
	 * @param  lin      Line
	 * @return          Value at _M[col + _dim[COL]*lin]
	 */
    T
    inline operator()          (const size_t col, const size_t lin) const {
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
    operator()           (const size_t col, const size_t lin) {
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
    operator()           (const size_t col, const size_t lin, const size_t slc) const {
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
    operator()           (const size_t col, const size_t lin, const size_t slc) {
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
    operator()           (const size_t col, 
						  const size_t lin, 
						  const size_t cha,
						  const size_t set,
						  const size_t eco = 0,
						  const size_t phs = 0,
						  const size_t rep = 0,
						  const size_t seg = 0,
						  const size_t par = 0,
						  const size_t slc = 0,
						  const size_t ida = 0,
						  const size_t idb = 0,
						  const size_t idc = 0,
						  const size_t idd = 0,
						  const size_t ide = 0,
						  const size_t ave = 0) { 

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
    operator()           (const size_t col, 
						  const size_t lin, 
						  const size_t cha,
						  const size_t set,
						  const size_t eco = 0,
						  const size_t phs = 0,
						  const size_t rep = 0,
						  const size_t seg = 0,
						  const size_t par = 0,
						  const size_t slc = 0,
						  const size_t ida = 0,
						  const size_t idb = 0,
						  const size_t idc = 0,
						  const size_t idd = 0,
						  const size_t ide = 0,
						  const size_t ave = 0) const { 
		
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
    

	/**
	 * @brief           Elementwise addition with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m + t
	 */
	template <class S> inline friend Matrix<T>    
	operator+  (const S s, Matrix<T> &m) {
		return   m + T(s);
	}


	/**
	 * @brief           Elementwise subtraction from scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          -(m - s)
	 */
	template <class S> inline friend Matrix<T>
	operator-  (const S s, Matrix<T> &m) {
		return -(m - T(s));
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	template <class S> inline friend Matrix<T>    
	operator*  (const S& s, Matrix<T> &m) { 
		return   m * s;
	}


	/**
	 * @brief           Elementwise equality with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m == s
	 */
	inline friend Matrix<bool> 
	operator== (T s, Matrix<T> m) {
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
	operator>= (T s, Matrix<T> m) {
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
	operator<= (T s, Matrix<T> m) {
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
	operator!= (T s, Matrix<T> m) {
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
	operator>  (T s, Matrix<T> m) {
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
	operator<  (T s, Matrix<T> m) {
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
	operator&  (Matrix<bool>& mb, Matrix<T>& m) {
		return   m &  mb;
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
    Res                 (const size_t i)                                const {return _res[i];}
    
    
    /**
     * @brief           Rresolution a given dimension.
     *
	 * @param   i       Dimension
     * @return          Resolution
     */
    inline float&          
    Res                 (const size_t i)                                 {return _res[i];}
    
    
    /**
     * @brief           Get size a given dimension.
     *
	 * @param   i       Dimension
     * @return          Dimension
     */
    inline size_t          
    Dim                 (const size_t i)                                const {return _dim[i];}
    
    
    /**
     * @brief           Get reference to size a given dimension.
     *
     * @return          Number of rows.
     */
    inline size_t&          
    Dim                 (const size_t i)                                 {return _dim[i];}
    
    
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
    Dim                 (const int i)                                const {return _dim[i];}
    
    
    /**
     * @brief           Get reference to size a given dimension.
     *
     * @return          Number of rows.
     */
    inline size_t&          
    Dim                 (const int i)                                 {return _dim[i];}
    
    
    /**
     * @brief           Reset all dimensions to values in dim
	 *
	 * @param  dim      New dimensions
     */
    inline void         
    Dim                 (const size_t* dim)                                     const {
        for (size_t i = 0; i<INVALID_DIM; i++)
            _dim[i] = dim[i];
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
     * @brief           Resize to dims, reallocate and zero repository. 
	 *                  Needs to becalled after any resize operation on dimensios. 
	 *                  (i.e. M.Dim(2) = 10; Reset();)
     *
     * @param  dim      New dimensions
     */
    inline void         
    Resize              (const size_t& m, const size_t& n)                                      {

		_dim[0] = m;
		_dim[1] = n;

		for (size_t i = 2; i < INVALID_DIM; i++)
			_dim[i] = 1;

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
            _dim[i] = dim[i];

		_M.resize(Size());
		
		Zero();

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
     * @brief           Reset. i.e. reallocate and set all fields = T(0)
     */
    inline void         
    Reset               ()                                      {

		_M.resize(Size());
		
		Zero();


    }
    

    /**
     * @brief           Reset. i.e. Set all fields = T(0)
     */
    inline void         
    Zero               ()                                      {

		for (size_t i = 0; i < Size(); i++)
			_M[i] = (T) 0;

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
     * @brief           Matrix product. i.e. this * M.
     *
     * @param  M        The factor.
     */
    template<class S> Matrix<T>           
    operator->*         (Matrix<S>& M);
    
    
    /**
     * @brief           Elementwise substruction of two matrices
     *
     * @param  M        Matrix substruent.
     */
    template <class S> Matrix<T>           
    operator-           (const Matrix<S> &M);
    
    
    /**
     * @brief           Elementwise substruction all elements by a scalar
     *
     * @param  s        Scalar substruent.
     */
    template <class S> Matrix<T>           
    operator-           (const S s);
    
    
    /**
     * @brief           ELementwise substraction and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
	 * @return          Result
     */
    template <class S>  Matrix<T>           
    operator-=          (const Matrix<S> &M);
    
    
    /**
     * @brief           ELementwise substration with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator-=         (const S s);
    
    
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
     * @param  M        Matrix additive.
     */
    template <class S> Matrix<T>           
    operator+          (const Matrix<S> &M);
    
    
    /**
     * @brief           Elementwise addition iof all elements with a scalar
     *
     * @param  s        Scalar additive.
     */
    template <class S> Matrix<T>           
    operator+           (const S s);
    
    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator+=          (const Matrix<S> &M);
    
    
    /**
     * @brief           ELementwise addition with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
	 * @return          Result
     */
    template <class S > Matrix<T>           
    operator+=          (const S s);
    
    

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
	 * @param  M        The operand
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
	 * @return          Matrix of false where elements are equal s and true else.
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
     * @param  p        Power.
	 * @return          Result
     */
    Matrix<T>           
    operator^           (const float p);
    
    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Power.
	 * @return          Result
     */
    Matrix<T>           
    operator^=          (const float p);
    
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
	Random             ();    


    /**
     * @brief           Elementwise multiplication. i.e. this .* M.
     *
     * @param  M        Factor matrix.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator*          (const Matrix<S> &M);


    /**
     * @brief           Elementwise multiplication with a scalar. i.e. this * m.
	 *
	 * @param  s        Factor scalar
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator*          (const S& s);

    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator*=         (const Matrix<S> &M);
    
    
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
    operator/          (const Matrix<S> &M);

    
    /**
     * @brief           Elementwise division by scalar. i.e. this * m.
	 *
     * @param  s        The divisor.
	 * @return          Result
     */
    template <class S> Matrix<T>           
    operator/           (const S s);
    
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
    operator/=         (const S s);
    
    
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
     * @brief           Scalar product (complex: conjugate first vector) using <a href="http://www.netlib.org/blas/">BLAS</a> routines CDOTU and DDOT
     *
     * @param  M        Factor
     */
	T
    dotc (Matrix<T>& M);
    
	T
    dotu (Matrix<T>& M);
    
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
Matrix<T>::LinSpace (const T& start, const T& space, const T& end) {

	assert (space != T(0));

	Matrix<T> res;
	size_t n;

	n   = (size_t) ceil ((end - start) / space);
	res = Matrix<T>::Zeros (n,1);

	res[0] = start;

	for (int i = 1; i < n; i++)
		res[i] = res[i-1] + space;

	return res;

}



#include "Matrix_Constructors.hpp"
#include "Matrix_Operators.cpp"
#include "Matrix_ICE.cpp"
//#include "Matrix_FFT.hpp"

#endif // __MATRIX_H__
