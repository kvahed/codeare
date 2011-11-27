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

#include "cycle.h"            // FFTW cycle implementation

#include <complex>
#include <assert.h>

#include <iostream>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <vector>

#define ICE_SHRT_MAX 4095


#ifdef HAVE_MAT_H
#include "mat.h"
#endif

#include "Ptr.h"

/**
 * @brief raw data
 */
typedef std::complex<float> cplx;


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
    const size_t       
    Import              (const IceAs* ias);


    /**
     * @brief           Continue import from IceAs
     *                   
     * @param  ias      IceAs containing data
     * @param  pos      Import data starting at position pos of own repository
     * @return          Amount of data read
     */
    const size_t       
    Import              (const IceAs* ias, const size_t pos);


    /**
     * @brief           Import with MDH
     *                   
     * @param  ias      IceAs containing data
     * @param  mdh      Measurement data header      
     * @return          Amount of data read
     */
    const size_t       
    Import              (const IceAs* ias, sMDH* mdh);


    /**
     * @brief           Export data to ias
     *                   
     * @param  ias      IceAs for data export
     * @return          Amount of data exported
     */
    const size_t         
    Export              (IceAs* ias) const;


    /**
     * @brief           Partially export data to ias 
     * 
     * @param  ias      IceAs for data export
     * @param  pos      Export data starting at position pos of our repository
     */
    const size_t
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
    inline const T*            
    Data                ()  const {
        return &_M[0];
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
	 * @brief          Multiplication operator with different class
	 * 
	 * 
	 */
	
	template<class S> Matrix<T>
	operator*          (const Matrix<S>& M);

	
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
	friend Matrix<T>    
	operator+  (T s, Matrix<T> &m) {
		return   m + s;
	}


	/**
	 * @brief           Elementwise subtraction from scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          -(m - s)
	 */
	friend Matrix<T>    
	operator-  (T s, Matrix<T> &m) {
		return -(m - s);
	}


	/**
	 * @brief           Elementwise multiplication with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m * s
	 */
	friend Matrix<T>    
	operator*  (T s, Matrix<T> &m) { 
		return   m * s;
	}


	/**
	 * @brief           Elementwise multiplication of inverse with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          s * m.Inv()
	 */
	friend Matrix<T>    
	operator/  (T s, Matrix<T> &m) {
		return   s * m.Inv();
	}


	/**
	 * @brief           Elementwise equality with scalar (lhs)
	 *
	 * @param  s        Scalar lhs
	 * @param  m        Matrix rhs
	 * @return          m == s
	 */
	friend Matrix<bool> 
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
	friend Matrix<bool> 
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
	friend Matrix<bool> 
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
	friend Matrix<bool> 
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
	friend Matrix<bool> 
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
	friend Matrix<bool> 
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
	friend Matrix<T>    
	operator&  (Matrix<bool>& mb, Matrix<T>& m) {
		return   m &  mb;
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
    Row                  (const size_t r) const;
	
    
    /**
     * @brief            Operates only on inner 2D: Get a Row of data
     *  
     * @param  c         Column
     * @return           Copy data into new vector
     */
    Matrix <T>           
    Column               (const size_t c) const;
	
    /**
     * @brief            Operates only on inner 3D: Get a slice of data
     *  
     * @param  s         Slice
     * @return           Copy data into new matrix and return
     */
    Matrix <T>           
    Slice                (const size_t s) const;
	
    /**
     * @brief            Operates only on inner 3D: Get a slice of data
     *  
     * @param  v         Volume         
     * @return           Copy data into new matrix and return
     */
    Matrix <T>           
    Volume               (const size_t v) const;
	

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
     * @brief           Get number of rows, i.e. tmp = size(this); tmp(1).
     *
     * @return          Number of rows.
     */
    inline size_t                 
    m                   () const {return _dim[0];}


    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    inline size_t                 
    n                   () const {return _dim[1];}


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
    

    /**
     * @brief           Get the number of matrix cells, i.e. dim_0*...*dim_16.
     *
     * @return          Number of cells.
     */
    const size_t
    Size                ()                                    const;
    
    
    /**
     * @brief           Check if we are XD (i.e. X dimensions > 1)
     *
	 * @param  d        Dimensions
     * @return          XD matrix?
     */
    const bool                
    IsXD                (const size_t d)                         const;
    
    
    /**
     * @brief           Check if we are 2D (i.e. COL, LIN)
     *
     * @return          2D matrix?
     */
    const bool                
    Is1D                ()                                    const;
    
    
    /**
     * @brief           Check if we are 2D (i.e. COL, LIN)
     *
     * @return          2D matrix?
     */
    const bool                
    Is2D                ()                                    const;
    
    
    /**
     * @brief           Check if we are 2D (i.e. COL, LIN)
     *
     * @return          2D matrix?
     */
    const bool                
    IsZero              ()                                    const;
    
    
    /**
     * @brief           Check if we are 3D (i.e. COL, LIN, SLC)
     *
     * @return          3D matrix?
     */
    const bool                
    Is3D                ()                                    const;
    
    
    /**
     * @brief           Check if we are 4D (i.e. COL, LIN, SLC)
     *
     * @return          3D matrix?
     */
    const bool 
    Is4D                ()                                    const;
    

	/**
	 * @brief           Check if empty
	 *
	 * @return          Empty: true;
	 */
	const bool 
	Empty               ()                                    const ;

    
    /**
     * @brief           Get the number of matrix cells, i.e. Size * sizeof(T).
     *
     * @return          Size in RAM in bytes.
     */
    const size_t                  
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
    Matrix<T>           
    operator-           (Matrix<T> &M);
    
    
    /**
     * @brief           Elementwise substruction all elements by a scalar
     *
     * @param  s        Scalar substruent.
     */
    Matrix<T>           
    operator-           (const T s);
    
    
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
	Random               ();    


    /**
     * @brief           Elementwise multiplication. i.e. this .* M.
     *
     * @param  M        Factor matrix.
	 * @return          Result
     */
    Matrix<T>           
    operator*           (const Matrix<T> &M);


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
	 * @return          Result
     */
    Matrix<T>           
    operator*=           (const Matrix<T> &M);
    
    
    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this *= s.
     *
     * @param  s        Factor scalar.
	 * @return          Result
     */
    Matrix<T>           
    operator*=           (const T s);
    
    
    /**
     * @brief           ELementwise division and assignment operator. i.e. this = this ./ M.
     *
     * @param  M        Divisor matrix.
	 * @return          Result
     */
    Matrix<T>           
    operator/=           (const Matrix<T> &M);
    
    
    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Divisor scalar.
	 * @return          Result
     */
    Matrix<T>           
    operator/=           (T s);
    
    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
	 * @return          Result
     */
    Matrix<T>           
    operator+=           (const Matrix<T> &M);
    
    
    /**
     * @brief           ELementwise addition with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
	 * @return          Result
     */
    Matrix<T>           
    operator+=           (T s);
    
    
    /**
     * @brief           ELementwise substraction and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
	 * @return          Result
     */
    Matrix<T>           
    operator-=           (const Matrix<T> &M);
    
    
    /**
     * @brief           ELementwise substration with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
	 * @return          Result
     */
    Matrix<T>           
    operator-=           (const T s);
    
    
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
     * @brief           Print dimensions to STL string.
     *
     * @return          Dimension string
     */
    const std::string       
    DimsToString        () const;
    

    /**
     * @brief           Print dimensions to C string.
     *
     * @return          Dimension string
     */
    const char*
    DimsToCString        () const;
    

    /**
     * @brief           Print resolutions to STL string.
     *
     * @return          Dimension string
     */
    const std::string       
    ResToString        () const;
    

    /**
     * @brief           Print resolutionss to C string.
     *
     * @return          Dimension string
     */
    const char*
    ResToCString        () const;
    

    /**
     * @brief           Primitive dump column-major to file.
     * 
     * @param  fname    File name.
     * @return          Success.
     */
    bool                
    PRDump              (const std::string fname) const;
    

    /**
     * @brief           Dump to <a href="http://www.hdfgroup.org/HDF5/" target="io">HDF5</a> file.
     * 
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    H5Dump              (const std::string fname, const std::string dname = "", const std::string dloc = "/") const;
    

    /**
     * @brief           Read from <a href="http://www.hdfgroup.org/HDF5/" target="io">HDF5</a> file.
     *
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    H5Read              (const std::string fname, const std::string dname = "", const std::string dloc = "/");
    

    /**
     * @brief           Dump to <a href="http://cdf.gsfc.nasa.gov/" target="io">CDF</a> file.
     * 
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    CDFDump             (const std::string fname, const std::string dname = "", const std::string dloc = "/") const;
    

    /**
     * @brief           Read from <a href="http://cdf.gsfc.nasa.gov/" target="io">CDF</a> file.
     *
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    CDFRead             (const std::string fname, const std::string dname = "", const std::string dloc = "/");
    

    /**
     * @brief           Dump to <a href="http://nifti.nimh.nih.gov/">NIFTI</a> file.
     * 
     * @param  fname    File name.
     * @return          Success.
     */
    bool                
    NIDump              (const std::string fname) const;
    

    /**
     * @brief           Read from <a href="http://nifti.nimh.nih.gov/">NIFTI</a> file.
     *
     * @param  fname    File name.
     * @return          Success.
     */
    bool                
    NIRead              (const std::string fname);
    

#ifdef HAVE_MAT_H

    /**
     * @brief           Dump to <a href="http://www.mathworks.com">MATLAB</a> file.
     * 
     * @param  file     File handle.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    MXDump              (MATFile* file, const std::string dname = "", const std::string dloc = "/") const;
    

    /**
     * @brief           Dump to <a href="http://www.mathworks.com">MATLAB</a> file.
     * 
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    MXDump              (const std::string fname, const std::string dname = "", const std::string dloc = "/") const;
    

    /**
     * @brief           Read from <a href="http://www.mathworks.com">MATLAB</a> file.
     *
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
     * @return          Success.
     */
    bool                
    MXRead              (const std::string fname, const std::string dname = "", const std::string dloc = "/");

#endif    

    /**
     * @brief           Dump matrix to file.
     * 
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
	 * @param  ios      IO strategy (HDF5, MATLAB, NIFTI, primitive)
     * @return          Success.
     */
    bool                
    Dump                (const std::string fname, const std::string dname = "", const std::string dloc = "/", const io_strategy ios = HDF5) const ;
    

    /**
     * @brief           Read matrix from file.
     *
     * @param  fname    File name.
	 * @param  dname    Dataset name.
	 * @param  dloc     Dataset location.
	 * @param  ios      IO strategy (HDF5, MATLAB, NIFTI, SYNGO)
     * @return          Success.
     */
    bool                
    Read                (const std::string& fname, const std::string& dname = "", const std::string& dloc = "/", const io_strategy& ios = HDF5);

    
	/**
	 * @brief           Read from <a href="http://www.medical.siemens.com/">Syngo MR</a> meas file
	 *
	 * @param  fname    File name
	 * @param  version  File version
	 * @return          Success
	 */
	bool
	RAWRead             (const std::string& fname, const std::string& version);

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
    Max                 ();    
    
    
    /**
     * @brief           Smallest element of the matrix. i.e. min(M0).
     *
     * @return          The scalar value of the smallest element.
     */
    T                   
    Min                 ();
    
    
	/**
	 * @brief           Conjugation w/o transposing
	 *
	 * @return          Conjugated
	 */
	Matrix<T>
	Conj                () const ;


	/**
	 * @brief           In-place conjugation w/o transposing
	 */
	void
	Conj                ();


	/**
	 * @brief           In-place Isotropic resampling up to 3D to
	 *
	 * @param  f        Resampling factor (i.e. 0.5 = each dimension halved)
	 * @param  i        Interpolation method @see InterpMethod
	 */
	void
	Resample            (const Matrix<double>& f, const InterpMethod& i);


	/**
	 * @brief           Isotropic resampling up to 3D to
	 *
	 * @param  f        Resampling factor (i.e. 0.5 = each dimension halved)
	 * @param  i        Interpolation method @see InterpMethod
	 */
	Matrix<T>
	Resample            (const Matrix<double>& f, const InterpMethod& i) const;


	/**
	 * @brief           Resampling up to 3D to
	 *
	 * @param  f        Resampling factors <1 reduction >1 expansion
	 * @param  im       Interpolation method @see InterpMethod
	 */
	void
	Resample            (const double& f, const InterpMethod& im);


	/**
	 * @brief           Resampling up to 3D to
	 *
	 * @param  f        Resampling factors <1 reduction >1 expansion
	 * @param  im       Interpolation method @see InterpMethod
	 */
	Matrix<T>
	Resample            (const double& f, const InterpMethod& im) const;


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
	 * @brief           Absolute values matrix
	 *
	 * @return          Absolute values
	 */
	Matrix<T>
	Abs                 () const {
		
		Matrix<T> res(_dim);
		
#pragma omp parallel default (shared) 
		{
			
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
			
			for (size_t i = 0; i < Size(); i++)
				res[i] = abs(_M[i]);

		}		

		return res;
		
	}
	

	/**
	 * @brief           Arguments matrix
	 *
	 * @return          Arguments
	 */
    template<class S> Matrix<S>           
	Arg                 () const;
		
	/**
	 * @brief           Real values matrix
	 *
	 * @return          Real values
	 */
    template<class S> Matrix<S>           
	Real                () const;
	

	/**
	 * @brief           Imaginary values matrix
	 *
	 * @return          Imaginary values
	 */
    template<class S> Matrix<S>           
	Imag                () const;
	

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
     * @brief           Fourier transform (<a href="fftw.org" target="fftw">FFTW</a>) over all dimensions.
     *
     * @return          Fourier transform.
     */
    Matrix<T>           
    FFT                 ()            const;
    
    /**
     * @brief           Inverse Fourier transform <a href="fftw.org" target="fftw">FFTW</a> over all dimensions
     *
     * @return          Inverse Fourier transform.
     */
    Matrix<T>           
    IFFT                ()            const;
    
    /**
     * @brief           MATLAB-like fftshift.<br/> 
	 *                  i.e. shift zero component to centre of volume for viewing
     *
     * @return          FFT shift.
     */
    Matrix<T>           
    FFTShift            (const size_t d = 0)            const;
    
    /**
     * @brief           MATLAB-like ifftshift.<br/>
	 *                  i.e. shift zero component to centre of volume for viewing. WORKS ONLY ON EVEN 
     *
     * @return          FFT shift.
     */
    Matrix<T>           
    IFFTShift           (const size_t d = 0)            const;
    

    /**
     * @brief           ND Hann window
     *
     * @return          FFT shift.
     */
    Matrix<T>           
    HannWindow          (const size_t d = 0)            const;
    

    /**
     * @brief           Sum of squares. 
     *
	 * @param  d        Dimension to eliminate. <br/>If not given, it is done over the outermost
     * @return          Sum of squares
     */
	Matrix<T>
    SOS                 (const size_t d = 0)           const;
    

    /**
     * @brief           Sum of squares on itself. 
     *
	 * @param  d        Dimension to eliminate. <br/>If not given, it is done over the outermost
     * @return          Sum of squares
     */
	void
    SOS                 (const size_t d = 0);
    

    /**
     * @brief           Mean along any dimension 
     *
	 * @param  d        Dimension to eliminate. If not given, done over the outermost
     * @return          Mean
     */
	Matrix<T>
    Mean                (const size_t d = 0)           const;
    

    /**
     * @brief           Mean of this along any dimension  
     *
	 * @param  d        Dimension to eliminate. If not given, done over the outermost
     * @return          Mean
     */
	void
    Mean                (const size_t d = 0);
    

    /**
     * @brief           Mean along any dimension 
     *
	 * @param  d        Dimension to eliminate. If not given, done over the outermost
     * @return          Mean
     */
	Matrix<T>
    Sum                 (const size_t d = 0)           const;
    

    /**
     * @brief           Mean of this along any dimension  
     *
	 * @param  d        Dimension to eliminate. If not given, done over the outermost
     * @return          Mean
     */
	void
    Sum                 (const size_t d = 0);
    

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
	 * @brief           Subscript of 1st dimension from index
	 *
	 * @param  ind      Index
	 * @return          Subscript of 1st dimension
	 */ 
	const size_t 
	Ind2i               (const size_t& ind) const;


	/**
	 * @brief           Subscript of 2nd dimension from index
	 *
	 * @param  ind      Index
	 * @return          Subscript of 2nd dimension
	 */ 
	const size_t 
	Ind2j               (const size_t& ind) const;


	/**
	 * @brief           Subscript of 3rd dimension from index
	 *
	 * @param  ind      Index
	 * @return          Subscript of 3rd dimension
	 */ 
	const size_t 
	Ind2k               (const size_t& ind) const;


	/**
	 * @brief           Subscript of 4th dimension from index
	 *
	 * @param  ind      Index
	 * @return          Subscript of 4th dimension
	 */ 
	const size_t 
	Ind2l               (const size_t& ind) const;


	/**
	 * @brief           Subscript of x-th dimension from index
	 *
	 * @param  ind      Index
	 * @param  dim      Dimension
	 * @return          Subscipt of x-th dimension
	 */ 
	const size_t 
	Ind2x               (const size_t& ind, const size_t& dim) const;


	/**
	 * @brief           Highest occupied dimension
	 *
	 * @return 
	 */
	size_t
	HDim                () const;


	/**
	 * @brief           Print dimension to std out
	 *
	 * @return
	 */
	void
	PrintDims           () const;


    /**
     * @brief           Inversion of positive definite matrix through LU factorisation.<br/>
	 *                  Wrapper to <a href="http://www.netlib.org/lapack/">LAPACK</a> drivers (xGETRI & xGERTF).
     *
     * @return          Inverse
     */
	Matrix<T>
    Inv   ()            const;
    

    /**
     * @brief           Moore-Penrose general inverse.<br/>
	 *                  Wrapper to <a href="http://www.netlib.org/lapack/">LAPACK</a> routines (xGELSD).
     *
     * @return          Pseudo-inverse
     */
	Matrix<T> 
    Pinv   ();
    

    /**
     * @brief           Cholesky factorisation of a positive definite matrix.<br/>
	 *                  Wrapper to <a href="http://www.netlib.org/lapack/">LAPACK</a> routines (xPOTRF).
	 * 
     * @param  uplo     Store upper or lower triangular matrix {'U','L'}
     * @return          Factorisation
     */
	Matrix<T> 
    Cholesky (const char uplo);
    

    /**
     * @brief           Euclidean norm using <a href="http://www.netlib.org/blas/">BLAS</a> routines xNRM2.
     *
     * @return          Norm
     */
	T
    Norm ()             const;
    

    /**
     * @brief           Scalar product (complex: conjugate first vector) using <a href="http://www.netlib.org/blas/">BLAS</a> routines CDOTU and DDOT
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
     * @brief           Compute eigen values with <a href="http://www.netlib.org/lapack/">LAPACK</a> routines xGEEV.
     *
	 * @param  ev       Vector containing the computed eigenvalues
	 * @param  cv       Compute also left hand eigen vectors
	 * @param  lev      Left hand eigen vectors 
	 * @param  rev      right hand eigen vectors 
     * @return          Feedback from Lapack operation
     */
    inline int
	EIG                 (Matrix<T>* ev, Matrix<T>* lev, Matrix<T>* rev, const bool cv = true);
    

    /**
     * @brief           Compute singular value decomposition with <a href="http://www.netlib.org/lapack/">LAPACK</a> routines xGESDD.
     *
	 * @param  lsv      Left hand singular vectors.
	 * @param  rsv      Right hand singular vectors.
	 * @param  sv       Sorted singular values.
	 * @param  jobz     @see http://www.netlib.org/lapack/double/dgesdd.f @see http://www.netlib.org/lapack/double/cgesdd.f
     * @return          Info from Lapack operation.
     */
    inline int
	SVD                 (Matrix<T>* lsv, Matrix<T>* rsv, Matrix<T>* sv, const char jobz = 'A');
    
	
    //@}
    
    

private:
    
    size_t              _dim[INVALID_DIM]; /// Dimensions
    float               _res[INVALID_DIM]; /// Resolutions
	std::vector<T>      _M;
	std::string         _name; 

    /**
     * @brief           Matrix matrix product with BLAS.
     *
	 * @param  M        Multiplie with
	 * @param  transb   Transpose M
	 *
     * @return          Product of this and M.
     */
    Matrix<T>           
    GEMM                (Matrix<T> &M, const char transb = 'N');


    /**
     * @brief           Matrix vector product with BLAS.
     *
	 * @param  M        Multiplie with
	 * @param  transa   Transpose M
	 *
     * @return          Product of this and M.
     */
	Matrix<T>           
    GEMV                (Matrix<T> &M, const char transa = 'N');


	/**
	 * @brief           Adjust and resize for Syngo read
	 *
	 * @param  fname    Syngo MR meas file name
	 * @return          Success
	 */
	bool
	RSAdjust            (const std::string& fname);
    
};


template <class T> inline Matrix<T> 
Matrix<T>::Volume (const size_t s) const {
    
	assert (Is4D());
    
    Matrix<T> res;

	for (size_t j = 0; j < 3; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	size_t nc = _dim[0]*_dim[1]*_dim[2];

	memcpy (&res[0], &_M[s * nc], nc * sizeof(T));

	return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::Slice (const size_t s) const {
    
	assert (Is3D());
    
    Matrix<T> res;

	for (size_t j = 0; j < 2; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	size_t nc = _dim[0]*_dim[1];

	memcpy (&res[0], &_M[s * nc], nc*sizeof(T));

	return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::Row (const size_t r)  const {

	assert (Is2D());
    
    Matrix<T> res;

	res.Dim(0) = _dim[1];
	res.Reset();

	for (size_t i = 0; i < _dim[1]; i++)
		res[i] = _M[r + i*_dim[0]];

	return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::Column (const size_t c) const {
    
    Matrix<T> res;

	res.Dim(0) = _dim[0];
	res.Reset();

	memcpy (&res[0], _M[c*_dim[0]], _dim[0] * sizeof(T));

	return res;

}


template <class T> inline const size_t
Matrix<T>::Size() const {
    
    long size = 1;
    
    for (size_t i = 0; i < INVALID_DIM; i++)
        size *= _dim[i];
    
    return size;
    
}


template <class T> inline const size_t 
Matrix<T>::SizeInRAM() const {
    
    return Size() * sizeof(T);
    
}


template <> inline short 
Matrix<short>::Max() {
	
    short max = _M[0];
	
    for (size_t i = 0; i < Size(); i++)
        if (_M[i] > max)
            max = _M[i];
	
    return max;
	
}


template <> inline cplx
Matrix<cplx>::Max() {
		
	cplx   max = cplx(0.0,0.0);
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
Matrix<T>::Maxabs() {

    T max = abs(_M[0]);

    for (size_t i = 0; i < Size(); i++)
        if (abs(_M[i]) > abs(max))
            max = abs(_M[i]);

    return abs(max);

}


template <class T> inline T
Matrix<T>::Min() {

    T min = _M[0];

    for (size_t i = 0; i < Size(); i++)
        if (_M[i] < min)
            min = _M[i];

    return min;

}


template <class T> inline T  
Matrix<T>::Minabs() {

    T old = fabs(_M[0]);

    for (size_t i = 0; i < Size(); i++)
        if (fabs(_M[i]) < old)
            old = fabs(_M[i]);

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
			if (typeid (T) == typeid (double))
				res.At(i,j) =           At(j,i);  // Transpose
			else
				res.At(i,j) = std::conj(At(j,i)); // Conjugate transpose

    return res;

}

template<> inline void 
Matrix<cplx>::Random () {
	
	srand (time(NULL));

	for (size_t i = 0; i < Size(); i++)
		_M[i] = cplx ((float) rand() / (float) RAND_MAX*2-1, (float) rand() / (float) RAND_MAX*2-1);
	
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

#include "Matrix_Constructors.hpp"
#include "Matrix_IO.cpp"
#include "Matrix_Lapack.cpp"
#include "Matrix_Operators.cpp"
#include "Matrix_BLAS.cpp"
#include "Matrix_ICE.cpp"
#include "Matrix_FFT.hpp"
#include "Matrix_Algorithms.hpp"
#include "Matrix_ITK.cpp"

#endif // __MATRIX_H__
