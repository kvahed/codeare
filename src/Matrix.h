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

#include <complex>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
/**
 * @brief raw data
 */
typedef std::complex<float> raw;

extern "C" {
	void cgeev_(char* jobvl, char* jobvr, int* n, raw*    a, int* lda, raw*    w,              raw*    vl, int* ldvl, raw*    vr, int* ldvr, raw*    work, int* lwork, float* rwork, int* info);
	void dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda, double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info);
}

/**
 * Short test if the matrix is a vector.
 */
# define VECT(M) assert((M)->width() == 1 || (M)->height() == 1);

/**
 * Return absolute value.
 */
# define ABS(A) (A > 0 ? A : -A)

/**
 * Return rounded value as in MATLAB.
 * i.e. ROUND( 0.1) ->  0
 *      ROUND(-0.5) -> -1
 *      ROUND( 0.5) ->  1
 */
# define ROUND(A) ( floor(A) + ((A - floor(A) >= 0.5) ? (A>0 ? 1 : 0) : 0))

/**
 * Defined for memory allocation testing.
 */
static int nb_alloc = 0;


/**
 * Old friend pi.
 */
# define PI  3.1415926535897931159979634685441851615906





/**
 * @brief   Matrix template
 *          This class intends to offer a simple interface for handling
 *          MR data in a simple way. Version 0.3 copes with SIEMENS MRIR.
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
     * @brief           Construct a new 1^16 matrix initialized with scalar s.
     *
     * @param  s        Scalar Value.
     */
    inline              
    Matrix              (T s) ;
    
    
    /**
     * @brief           Constructs a 1^16 matrix matrix with desired dimensions values = 0.
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
    Matrix               (int col, int lin, int cha, int set, 
                          int eco, int phs, int rep, int seg, 
                          int par, int slc, int ida, int idb, 
                          int idc, int idd, int ide, int ave);
    
    
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

    // If compiled within IDEA we know of access specifiers.
	#ifdef PARC_MODULE_NAME
    
    /**
     * @brief           Reset and fill data from IceAs
     *                   
     * @param  ias      IceAs containing data
     * 
     * @return          Amount of data read
     */
    inline long         
    Import              (IceAs ias);

    /**
     * @brief           Import an amount of data from IceAs
     *                   
     * @param  ias      IceAs containing data
     * @param  pos      Import data starting at position pos own repository
     * 
     * @return          Amount of data read
     */
    inline long         
    Import              (IceAs ias, long pos);

    /**
     * @brief           Export data to ias
     *                   
     * @param  ias      IceAs for data export
     * 
     * @return          Amount of data exported
     */
    inline long         
    Export              (IceAs ias);

    /**
     * @brief           Partially export data to ias 
     * 
     * @param  ias      IceAs for data export
     * @param  pos      Export data starting at position pos of our repository
     */
    inline long         
    Export              (IceAs ias, long pos);
 
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
    operator[]          (int p)                             const;
    
    
    /**
     * @brief           Reference to pth element from data repository.
     *
     * @param  p        Requested position.
     * @return          Reference to _M[p].
     */
    T                   
    &operator[]         (int p)                              ;

    
    /**
     * @brief           Get element at position 
     *  
     * @param  pos      Position
     * @param  lin      Line
     *
     * @return          Value at _M[pos]
     */
    inline T            
    at                  (int pos)  const {
        return _M[pos];
    }

    
    /**
     * @brief            Reference to value at position
     *  
     * @param  pos       Position
     * @param  lin       Line
     *
     * @return           Reference to _M[pos]
     */
    inline T&           
    at                  (int pos) {
        return _M[pos];
    }

    
    /**
     * @brief           Get value in slice
     *  
     * @param  col      Column
     * @param  lin      Line
     *
     * @return          Value at _M[col + _dim[LIN]*lin]
     */
    inline T            
    at                  (int col, int lin)  const {
        return _M[col + _dim[LIN]*lin ];
    }

    
    /**
     * @brief            Reference to value in slice
     *  
     * @param  col       Column
     * @param  lin       Line
     *
     * @return           Reference to _M[col + _dim[LIN]*lin]
     */
    inline T&           
    at                  (int col, int lin) {
        return _M[col + _dim[LIN]*lin ];
    }

    
    /**
     * @brief            Get value in volume
     *  
     * @param  col       Column
     * @param  lin       Line
     * @param  slc       Slice
     *
     * @return           Value at _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc]
     */
    inline T            
    at                   (int col, int lin, int slc)  const {
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
    at                   (int col, int lin, int slc) {
        return _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc];
    }
    
    
    /**
     * @brief            Get value in volume
     *  
     * @param  col       Column
     * @param  lin       Line
     * @param  slc       Slice
     *
     * @return           Value at _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc]
     */
    inline T            
    at                   (int col, 
						  int lin, 
						  int cha,
						  int set,
						  int eco,
						  int phs,
						  int rep,
						  int seg,
						  int par,
						  int slc,
						  int ida,
						  int idb,
						  int idc,
						  int idd,
						  int ide,
						  int ave) const {
        return _M[
				  col +
				  _dim[COL]*lin +
				  _dim[COL]*_dim[LIN]*lin +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA] +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO] +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida*_dim[IDB]*idb +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida*_dim[IDB]*idb*_dim[IDC]*idc +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida*_dim[IDB]*idb*_dim[IDC]*idc*_dim[IDD]*idd +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida*_dim[IDB]*idb*_dim[IDC]*idc*_dim[IDD]*idd*_dim[IDE]*ide +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida*_dim[IDB]*idb*_dim[IDC]*idc*_dim[IDD]*idd*_dim[IDE]*ide*_dim[AVE]*ave];
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
    at                   (int col, 
						  int lin, 
						  int cha,
						  int set,
						  int eco,
						  int phs,
						  int rep,
						  int seg,
						  int par,
						  int slc,
						  int ida,
						  int idb,
						  int idc,
						  int idd,
						  int ide,
						  int ave) {
        return _M[
				  col +
				  _dim[COL]*lin +
				  _dim[COL]*_dim[LIN]*lin +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA] +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO] +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida*_dim[IDB]*idb +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida*_dim[IDB]*idb*_dim[IDC]*idc +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida*_dim[IDB]*idb*_dim[IDC]*idc*_dim[IDD]*idd +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida*_dim[IDB]*idb*_dim[IDC]*idc*_dim[IDD]*idd*_dim[IDE]*ide +
				  _dim[COL]*_dim[LIN]*lin*_dim[CHA]*_dim[ECO]*_dim[PHS]*phs*_dim[REP]*rep*_dim[SEG]*seg*_dim[PAR]*par*_dim[SLC]*slc*_dim[IDA]*ida*_dim[IDB]*idb*_dim[IDC]*idc*_dim[IDD]*idd*_dim[IDE]*ide*_dim[AVE]*ave];
    }
    
     /**
     * @brief           Get the element at position p of the vector, i.e. this(p).
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    T                  
    operator()          (int p) const;

    
    /**
     * @brief           Get value of pth element of repository.
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    T&                 
    operator()          (int p) ;

    
    /**
	 * @brief           Get value in slice
	 *
	 * @param  col      Column
	 * @param  lin      Line
	 *
	 * @return          Value at _M[col + _dim[LIN]*lin]
	 */
    T
    operator()          (int col, int lin) const {
        return _M[col + _dim[LIN]*lin ];
    }
    

    /**
	 * @brief           Reference to value in slice
	 *
	 * @param  col      Column
	 * @param  lin      Line
	 *
	 * @return          Reference to _M[col + _dim[LIN]*lin]
	 */
    T&                  
    operator()           (int col, int lin) {
        return _M[col + _dim[LIN]*lin ];
    }
    
    
    /**
     * @brief            Get value in volume
     *
     * @param  col       Column
     * @param  lin       Line
     * @param  slc       Slice
     *
     * @return           Value at _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc]
     */
    T                  
    operator()           (int col, int lin, int slc) const {
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
    operator()           (int col, int lin, int slc) {
           return _M[col + _dim[COL]*lin + _dim[COL]*_dim[LIN]*slc];
    }
    
    //@}
    




    /**
     * @name            Partial access functions.
     *                  Functions for access to parts of data.
     */
    
    //@{
    
    /**
     * @brief           Get a vector with as copy of line l of this Matrix, i.e. this(l,:).
     *                  NOT IMPLEMETED YET
     *
     * @param  lin      Requested line.
     * @return          Vector copied from requested row.
     */
    Matrix<T>          
    lin                 (int lin)                             const;
    
    
    /**
     * @brief           Get a matrix as a copy of lines l[i] of this matrix, i.e. this(l,:);
     *                  NOT IMPLEMETED YET
     *
     * @param  lin      Requested rows.
     * @return          Matrix copied from requested rows.
     */
    Matrix<T>           
    lin                 (Matrix<int> lin)                     const;
    
    
    /**
     * @brief           Get a vector with as copy of column c of this Matrix, i.e. this(:,c).
	 *                  NOT IMPLEMETED YET
     *
     * @param  col      Requested column.
     * @return          Vector copied from requested column.
     */
    Matrix<T>           
    col                 (int col)                             const;
    
    
    /**
     * @brief           Get a matrix as a copy of columns c[i] of this matrix, i.e. this(:,c);
     *                  NOT IMPLEMETED YET
     *
     * @param  col      Requested rows.
     * @return          Matrix copied from requested columns.
     */
    Matrix<T>           
    col                 (Matrix<int> col)                     const;


    /**
     * @brief           Get a vector with as copy of line l of this Matrix, i.e. this(l,:).
     *                  NOT IMPLEMETED YET
     *
     * @param  lin      Requested line.
     * @return          Vector copied from requested row.
     */
    Matrix<T>          
    slc                 (int slc)                             const;
    
    
    /**
     * @brief           Get a matrix as a copy of lines l[i] of this matrix, i.e. this(l,:);
     *                  NOT IMPLEMETED YET
     *
     * @param  lin      Requested rows.
     * @return          Matrix copied from requested rows.
     */
    Matrix<T>           
    slc                 (Matrix<int> slc)                     const;
    
    
    /**
	 * @brief           Get one or more elements of a vector, i.e. this(p).
	 *                  NOT IMPLEMENTED YET
	 *
	 * @param  p        Requested positions
	 * @return          Vector of requested elements.
	 */
	Matrix<T>
	operator()          (Matrix<int> &p)                    const;


	/**
	 * @brief           Get one or more elements of a matri, i.e. this=[1 2 3;4 5 6;7 8 9],
	 *                  m1=[2], m2=[1 3], returns [2 8]. NOT IMPLEMENTED YET.
	 *
	 * @param  col      Requested columns.
	 * @param  lin      Requested lines.
	 * @return          matrix of requested elements.
	 */
	Matrix<T>
	operator()          (Matrix<int> &col, Matrix<int> &lin)  const;




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
    int                 
    height              () const {return _dim[LIN];}
    
    
    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    int                 
    width               () const {return _dim[COL];}
    
    
    /**
     * @brief           Get number of rows, i.e. tmp = size(this); tmp(1).
     *
     * @return          Number of rows.
     */
    int                 
    m                   () const {return _dim[LIN];}


    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    int                 
    n                   () const {return _dim[COL];}


    /**
     * @brief           Get size a given dimension.
     *
     * @return          Number of rows.
     */
    inline int          
    Dim                 (int i)                                const {return _dim[i];}
    
    
    /**
     * @brief           Get size a given dimension.
     *
     * @return          Number of rows.
     */
    inline int&          
    Dim                 (int i)                                 {return _dim[i];}
    
    
    /**
     * @brief           Get size a given dimension.
     *
     * @return          Number of rows.
     */
    inline const int*   
    Dim                 ()                                     const {return _dim;}
    

    /**
     * @brief           Get size a given dimension.
     *
     * @return          Number of rows.
     */
    inline void         
    Dim                 (const int* dim)                                     const {
        for (int i = 0; i<INVALID_DIM; i++)
            _dim[i] = dim[i];
    }
    
    
    /**
     * @brief           Resize to dim and zero
     *
     * @param  dim      New dimensions
     */
    inline void         
    Reset               (const int* dim)                                      {

    	for (int i = 0; i < INVALID_DIM; i++)
            _dim[i] = dim[i];

        if (nb_alloc) {
            delete [] (_M);
            nb_alloc = 0;
        }

        _M = new T[Size()]();
        nb_alloc = 1;

    }
    

    /**
     * @brief           Resize to dim and zero
     *
     * @param  dim      New dimensions
     */
    inline void         
    Reset               ()                                      {

        if (nb_alloc) {
            delete [] (_M);
            nb_alloc = 0;
        }

        _M = new T[Size()]();

		for (int i = 0; i < Size(); i++)
			_M[i] = (T) 0;

        nb_alloc = 1;

    }
    

    /**
     * @brief           Get the number of matrix cells, i.e. dim_0*...*dim_16.
     *
     * @return          Number of cells.
     */
    long                
    Size                ()                                    const;
    
    
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
     * @brief           Substruct from 0. i.e. 0 - this.
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
     */
    Matrix<T>           
    operator+           ();
    
    
    /**
     * @brief           Transposition / Complex conjugation. i.e. this'.
     */
    Matrix<T>           
    operator!           () const;
    
    
    /**
     * @brief           Return a matrix with result[i] = (m[i] ? this[i] : 0).
     */
    Matrix<T>           
    operator&           (Matrix<bool> &M);
    
    
    /**
     * @brief           Scalar equality. result[i] = (this[i] == m).
     *
     * @param  s        Comparing scalar.
     * @return          True if and only if all elements of the matrix equal s.
     */
    Matrix<bool>        
    operator==          (T s);
    
    
    /**
     * @brief           Scalar inequality. result[i] = (this[i] != m). i.e. this ~= m
     *
     * @param  s        Comparing scalar.
     */
    Matrix<bool>        
    operator!=          (T s);
    
    
    /**
     * @brief           Scalar greater comaprison, result[i] = (this[i] > m). i.e. this > m
     *
     * @param  s        Comparing scalar.
     */
    Matrix<bool>        
    operator>           (T s);
    
    
    /**
     * @brief           Scalar greater or equal comparison. result[i] = (this[i] >= m). i.e. this >= m
     *
     * @param  s        Comparing scalar.
     */
    Matrix<bool>        
    operator>=          (T s);
    
    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] <= m). i.e. this <= m
     *
     * @param  s        Comparing scalar.
     */
    Matrix<bool>        
    operator<=          (T s);
    
    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] < m). i.e. this < m
     *
     * @param  s        Comparing scalar.
     */
    Matrix<bool>        
    operator<           (T s);
    
    
    /**
     * @brief           Elementwise equality, result[i] = (this[i] == m[i]). i.e. this == m
     *
     * @param  M        Comparing matrix.
     */
    Matrix<bool>        
    operator==          (Matrix<T> M);
    
    
    /**
     * @brief           Elementwise equality, result[i] = (this[i] != m[i]). i.e. this ~= m
     *
     * @param  M        Comparing matrix.
     */
    Matrix<bool>        
    operator!=          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] >= m[i]). i.e. this >= m
     *
     * @param  M        Comparing matrix.
     */
    Matrix<bool>        
    operator>=          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] <= m[i]). i.e. this <= m.
     *
     */
    Matrix<bool>        
    operator<=          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] > m[i]). i.e. this > m.
     *
     * @param  M        Comparing matrix.
     */
    Matrix<bool>        
    operator>           (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] < m[i]). i.e. this < m.
     *
     * @param  M        Comparing matrix.
     */
    Matrix<bool>        
    operator<           (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, tel que result[i] = (m[i] || this[i] ? 1 : 0). i.e. this | m.
     *
     * @param  M        Comparing matrix.
     */
    Matrix<T>           
    operator||          (Matrix<T> M);
    
    
    /**
     * @brief           Matrix comparison, tel que result[i] = (m[i] && this[i] ? 1 : 0). i.e. this & m.
     *
     * @param  M        Comparing matrix.
     */
    Matrix<T>           
    operator&&          (Matrix<T> &M);


    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Exponent.
     */
    Matrix<T>           
    operator^           (int p);
    
    //@}





    /**
     * @name Data manipulation functions.
     *       Data manipulation functions.
     */
    
    //@{
    
    /**
     * @brief           Fill with random data
     */
    /**
     * @brief           Resize to dim and zero
     */
    inline void         
	Random               ();    

    /**
     * @brief           Matrix multiplication. i.e. this * M.
     *
     * @param  M        The factor.
     */
    Matrix<T>           
    operator*           (Matrix<T> &M);
    
    /**
     * @brief           Elementwise multiplication with a scalar. i.e. this * m.
     */
    Matrix<T>           
    operator*           (T s);
    
    /**
     * @brief           Matrix multiplication. i.e. this * M.
     *
     * @param  M        The factor.
     */
    Matrix<T>           
    operator/           (Matrix<T> &M);
    
    /**
     * @brief           Elementwise multiplication with a scalar. i.e. this * m.
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
     * @brief           Dump binary matrix column-major.
     * 
     * @param  fname    File name.
     * @return          Success.
     */
    bool                
    dump                (const char* fname);
    

    /**
     * @brief           Read in binary matrix column-major.
     *
     * @param  fname    File name.
     * @return          Success.
     */
    bool                
    read                (const char* fname);
    
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
     * @brief           Transposition.
     *
     * @return          The transposed matrix
     */
    Matrix<T>           
    tr()                const;
    
    //@}


    /**
     * @name Lapack functions.
     *       Only available when lapack detected.
     */
    
    //@{
    

    /**
     * @brief           Compute eigen values.
     *
	 * @param  out      Vector containing the computed eigenvalues.
     * @return          Info from Lapack operation.
     */
    inline int
	geev                (Matrix<raw>& out);
    
    
    //@}
    
    

private:
    
    int                 _dim[INVALID_DIM]; /// Dimnesions
    T*                  _M;                /// Data repository
    
};


template <class T> 
Matrix<T>::Matrix () {

    for (int i = 0; i < INVALID_DIM; i++)
        _dim [i] = 1;

    _M = new T[1];
    nb_alloc = 1;


};


template <class T> 
Matrix<T>::Matrix (T s) {

    for (int i = 0; i < INVALID_DIM; i++)
        _dim [i] = 1;

    _M = new T[1];
    nb_alloc = 1;

    _M[0] = s;

};


template <class T>
Matrix<T>::Matrix (int col, int lin, int cha, int set, 
                   int eco, int phs, int rep, int seg, 
                   int par, int slc, int ida, int idb, 
                   int idc, int idd, int ide, int ave) {

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

    _M = new T[Size()]();
    nb_alloc = 1;


};


template <class T>
Matrix<T>::Matrix (const Matrix<T> &M)
{

	for (int i = 0; i < INVALID_DIM; i++) 
		_dim[i] = M.Dim(i);

	_M = new T[Size()]();

	for (int k = 0; k < M.Size(); k++)
		_M[k] = M[k];
	
	nb_alloc = 1;

}

#ifdef PARC_MODULE_NAME

template <class T> long 
Matrix<T>::Import     (IceAs ias, long pos) {
    
    ICE_SET_FN("Matrix<T>::Import(IceAs, long)")
        
    int  i    = 0;
    long size = 1;
    
    for (i = 0; i < INVALID_DIM; i++)
        size *= (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));
    
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (i = 0; i < size; i++, data++)
        _M[i+pos] = *data;
    
    return size;
    
}


template <class T> long 
Matrix<T>::Import(IceAs ias) {
    
    ICE_SET_FN("Matrix<T>::Import(IceAs)")
        
        int i;
    
    for (i = 0; i < INVALID_DIM; i++)
        _dim[i] = (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));
    
    _M = new T[Size()]();
    nb_alloc = 1;
    
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (i = 0; i < Size(); i++, data++)
        _M[i] = *data;
    
    return Size();
    
}


template <class T> long 
Matrix<T>::Export(IceAs ias) {
    
    ICE_SET_FN("Matrix<T>::Export(IceAs)")
        
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (int i = 0; i < Size(); i++, data++)
        *data = _M[i];
    
    return Size();
    
}


template <class T> long
Matrix<T>::Export(IceAs ias, long pos) {

    ICE_SET_FN("Matrix<T>::Export(IceAs, long)")
        
        int  i    = 0;
    long size = 1;
    
    for (i = 0; i < INVALID_DIM; i++) {
        size *= (ias.getLen(IceDim(i)) <= 1) ? 1 : ias.getLen(IceDim(i));
    }
    
    T* data = (T*) ias.calcSplObjStartAddr() ;
    
    for (i = 0; i < size; i++, data++)
        *data = _M[i+pos];
    
    return size;
    
}

#endif


template <class T> 
Matrix<T>::~Matrix() {
    
#ifdef PARC_MODULE_NAME
    ICE_SET_FN ("Matrix<T>::~Matrix()")
    ICE_WARN   ("Freeing " << (float)Size() * sizeof(T) / 1024 << " kB of RAM.");
#endif

    if (nb_alloc) {
    	delete [] (_M);
        nb_alloc = 0;
    }
    
}


template <class T> Matrix<T> 
Matrix<T>::operator*(Matrix<T> &M) {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	for (int i = 0; i < Size(); i++)
		res[i] = _M[i] * M[i];

	return res;

}


template <class T> Matrix<T> 
Matrix<T>::operator* (T s) {
    
    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	for (int i = 0; i < Size(); i++)
		res[i] = _M[i] * s;

	return res;

}


template <class T> Matrix<T> 
Matrix<T>::operator/(Matrix<T> &M) {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	for (int i = 0; i < Size(); i++)
		(M[i] != (T)0) ? res[i] = _M[i] / M[i] : 0;

	return res;

}


template <class T> Matrix<T> 
Matrix<T>::operator/ (T s) {
    
	assert (s != (T)0);

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	for (int i = 0; i < Size(); i++)
		res[i] = _M[i] * s;

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


template <class T> T           
Matrix<T>::operator[]  (int p) const {
    
    assert(p >= 0);
    assert(p <  Size());
    
    return _M[p];
    
}


template <class T> T           
&Matrix<T>::operator[] (int p) {
    
    assert(p >= 0);
    assert(p <  Size());
    
    return _M[p];
    
}


template <class T> Matrix<T>   
Matrix<T>::operator()(Matrix<int> &m) const {
    
    VECT(&m);
    
    Matrix<T> res(1, m.Size());
    
    for (int i = 0; i < res.Size(); i++)
        res[i] = (*this)[m[i]];
    
    return res;

};


template <class T> Matrix<T> 
Matrix<T>::operator() (Matrix<int> &m1, Matrix<int> &m2) const {

    Matrix<T> res(m1.Size(),  m2.Size());
    
    for (int i = 0; i < m1.Size(); i++)
        for (int j = 0; j < m2.Size(); j++)
            res[i * m2.Size() + j] = _M[m1[i] * _dim[LIN] + m2[j]];
    
    return res;
	
}


template <class T>
T Matrix<T>::operator() (int a) const {

    assert(a >= 0);
    assert(a <  Size());

    return _M[a];

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
		float tmp = 0.0;
		
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
Matrix<T> Matrix<T>::operator=(Matrix<T> &M) {
    
    int i;

    for (i = 0; i < INVALID_DIM; i++) 
        _dim[i] = M.Dim()[i];
    
    if (nb_alloc) {
        delete[](_M);
        nb_alloc = 0;
    }

    _M = new T[Size()]();
    nb_alloc = 1;

    for (i = 0; i < Size(); i++)
        _M[i] = M[i];

    return *this;

}


template <class T>
Matrix<bool> Matrix<T>::operator==(T s)    {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] == s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator>=(T s) {

    Matrix<bool> res(_dim);
    
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] >= s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator<=(T s) {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] <= s);
    
    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator!=(T s) {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] != s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator<(T s)    {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] < s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator>(T s) {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] > s);

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator->*(Matrix<T> &M) {

    return this->prod(M);

}

template <class T>
Matrix<T> Matrix<T>::operator!() const {

    return this->tr();

}

template <class T>
Matrix<T> Matrix<T>::operator&(Matrix<bool> &M)    {

    for (int i = 0; i < INVALID_DIM; i++) 
        assert (_dim[i] == (int) M.Dim()[i]);


    int k = 0, i;
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

template <class T>
Matrix<T> Matrix<T>::operator&&(Matrix<T> &M) {

    assert(M.Size() == Size());

    Matrix<T> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (M[i] && _M[i]);

    return res;

}

template <class T>
Matrix<T> Matrix<T>::operator||(Matrix<T> M) {

    assert(M.Size() == Size());

    Matrix<T> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (M[i] || _M[i]);

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator==(Matrix<T> M) {
    assert(_dim[COL] == M.width() && _dim[LIN] == M.height());
    Matrix<bool> res(_dim);
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] == M[i]);
    return res;
}

template <class T>
Matrix<bool> Matrix<T>::operator>=(Matrix<T> M) {
    assert(_dim[COL] == M.width() && _dim[LIN] == M.height());
    Matrix<bool> res(_dim);
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] >= M[i]);
    return res;
}

template <class T>
Matrix<bool> Matrix<T>::operator<=(Matrix<T> M) {
    assert(_dim[COL] == M.width() && _dim[LIN] == M.height());
    Matrix<bool> res(_dim);
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] >= M[i]);
    return res;
}

template <class T>
Matrix<bool> Matrix<T>::operator!=(Matrix<T> M) {
    assert(_dim[COL] == M.width() && _dim[LIN] == M.height());
    Matrix<bool> res(_dim);
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] != M[i]);
    return res;
}

template <class T>
Matrix<bool> Matrix<T>::operator>(Matrix<T> M) {
    assert(_dim[COL] == M.wridth() && _dim[LIN] == M.height());
    Matrix<bool> res(_dim);
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] > M[i]) ? true : false;
    return res;
}

template <class T>
Matrix<bool> Matrix<T>::operator<(Matrix<T> M) {

    assert(_dim[COL] == M.width() && _dim[LIN] == M.height());

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] < M[i]);

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator-() {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

    for (int i = 0; i < Size(); i++)
        res[i] =- _M[i];

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator-(Matrix<T> &M) {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	for (int i = 0; i < Size(); i++)
		res[i] = _M[i] - M[i];

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator-(T s) {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	for (int i = 0; i < Size(); i++)
		res[i] = _M[i] - s;

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator+() {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

    for (int i = 0; i < Size(); i++)
        res[i] =+ _M[i];

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator+(Matrix<T> &M) {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	for (int i = 0; i < Size(); i++)
		res[i] = _M[i] + M[i];

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator+(T s) {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

	for (int i = 0; i < Size(); i++)
		res[i] = _M[i] + s;

	return res;

}


template <class T>
static T power(T p, int b) {

    assert(b > 0);

    T res = p;

    b--;

    while (b > 0) {
        res *= p;
        b--;
    }

    return res;

}

template <class T>
Matrix<T> Matrix<T>::operator^(int p) {
    Matrix<T> res(_dim[LIN], _dim[COL]);
    for (int i = 0; i < Size(); i++)
        if (p == 0)
            res[i] = 1;
        else
            if (p > 0)
                res[i] = power(_M[i], p);
            else
                res[i] = 1 / power(_M[i], -p);
    return res;
}


template <class T>
Matrix<T> Matrix<T>::prod(Matrix<T> &m) {
    
    assert(_dim[COL] == m.height());
    Matrix<T> res(_dim[LIN], m.width());

    for (int i = 0; i < res.height(); i++)
        for (int j = 0; j < res.width(); j++) {
            res[i * res.width() + j] = 0;
            for (int k = 0; k < _dim[COL]; k++)
                    res[i * res.width() + j] += _M[i * _dim[COL] + k] *
                    m[k * m.width() + j];
        }
    return res;
}

template <class T>
Matrix<T> Matrix<T>::tr() const {

    Matrix<T> res(_dim);

    for (int i = 0; i < res.height(); i++)
        for (int j = 0; j < res.width(); j++)
            res[i * res.width() + j] = _M[j * _dim[COL] + i];

    return res;

}

template <class T>
bool Matrix<T>::dump (const char* fname) {

	// TODO: Error checking.
	//       File not found.
	//       HDF5 file format.

    if (fname != "") {

        std::ofstream fout(fname , std::ios::out | std::ios::binary);

        for (int i = 0; i < Size(); i++) {
			std::cout << _M[i] << std::endl;
            fout.write ((char*)(&(_M[i])), sizeof(T));
		}
		        
        fout.close();

    }

	return true;
    
}

template <class T>
bool Matrix<T>::read (const char* fname) {

	// TODO: Error checking.
	//       File not found.
	//       EOF before Size.
	//       HDF5 file format.

    if (fname != "") {

        std::ifstream fin  (fname , std::ios::in | std::ios::binary);

        for (int i = 0; i < Size(); i++)
            fin.read  ((char*)(&(_M[i])), sizeof(T));
        
        fin.close();

    }

	return true;
    
}


template<> inline std::ostream&  
Matrix<short>::Print     (std::ostream &os) const {

	int i,j;
	
	for (i = 0; i < _dim[COL]; i++) {
		for(j = 0; j < _dim[LIN]; j++)
			printf ("%4i ", _M [i * _dim[COL] + j]);
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<double>::Print     (std::ostream &os) const {

	int i,j;
	
	for (i = 0; i < _dim[COL]; i++) {
		for(j = 0; j < _dim[LIN]; j++)
			printf ("%.2f ", _M [i * _dim[COL] + j]);
		printf("\n");
	}
	
	return os;
	
}


template<> inline std::ostream&  
Matrix<raw>::Print       (std::ostream& os) const {
	
	int i,j;
	
	for (i = 0; i < _dim[COL]; i++) {
		for(j = 0; j < _dim[LIN]; j++)
			printf ("%.2f+%.2fi ", _M [i*_dim[COL]+j].real(), _M [i*_dim[COL]+j].imag());
		printf("\n");
	}
	
	return os;

}


template <class T> std::ostream& 
operator<<    (std::ostream& os, Matrix<T>& M) {

	M.Print(os);
	return os;

}


template <>
inline int Matrix<double>::geev (Matrix<raw>& out) {
	
	char    jobvl = 'N';
	char    jobvr = 'N';
	
	int     n     = _dim[COL];

	double* a     = new double[Size()];

	for (int i = 0; i < Size(); i++)
		a[i] = _M[i];

	int     lda   = n;
	
	double* wr    = new double[n];
	double* wi    = new double[n];
	double* vl    = new double[1];
	double* vr    = new double[1];

    int     ldvr  =  1;
    int     ldvl  =  1;
	
	double* work  = new double[1];

	int     info  =  0;
	int     lwork = -1;

	dgeev_ (&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
	lwork = (int) work[0];

	//std::cout << "cgeev_ (lwork = " << lwork << ")" << std::endl;

	delete [] work;
	work = new double[lwork];

	dgeev_ (&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);

	//std::cout << "dgeev_ (info = " << info << ")" << std::endl;
	
	out.Dim(COL) = n;
	out.Reset();
	
	for (int j = 0; j < n; j++)
		out[j] = raw (wr[j], wi[j]);

	std::cout << out << std::endl;


	delete [] a;
	delete [] wr;
	delete [] wi;
	delete [] vl;
	delete [] vr;
	delete [] work;

	return info;

}


template <>
inline int  Matrix<raw>::geev (Matrix<raw>& out) {
	
	char    jobvl = 'N';
	char    jobvr = 'N';
	
	int     n     = _dim[COL];

	raw*    a     = new raw[Size()];

	for (int i = 0; i < Size(); i++)
		a[i] = _M[i];

	int     lda   = n;
	
	raw*    w     = new raw[n];
	raw*    vl    = new raw[1];
	raw*    vr    = new raw[1];

    int     ldvr  =  1;
    int     ldvl  =  1;
	
	raw*    work  = new raw[1];

	int     info  =  0;
	int     lwork = -1;

	float*  rwork = new float[2*n];

	cgeev_ (&jobvl, &jobvr, &n, _M, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
	lwork = (int) work[0].real();

	//std::cout << "cgeev_ (lwork = " << lwork << ")" << std::endl;

	delete [] work;
	work = new raw[lwork];

	cgeev_ (&jobvl, &jobvr, &n, _M, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

	//std::cout << "cgeev_ (info = " << info << ")" << std::endl;

	out.Dim(COL) = n;
	out.Reset();
	
	for (int j = 0; j < n; j++)
		out[j] = w[j];

	std::cout << out << std::endl;

	delete [] w;
	delete [] vl;
	delete [] vr;
	delete [] work;

	return info;

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
		_M[i] = (short) 4096 * (double)rand() / (double)RAND_MAX*2-1;

}
    

#endif // __MATRIX_H__


