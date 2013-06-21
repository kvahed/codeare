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

#endif

#include "config.h"
#include "common.h"

#include "OMP.hpp"
#include "Complex.hpp"
#include "Container.hpp"

#include "cycle.h"            // FFTW cycle implementation

#include <assert.h>
#include <iostream>
#include <memory>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <ostream>
#include <string>
#include <cstring>
#include <algorithm>

#define ICE_SHRT_MAX 4095

#ifdef HAVE_MAT_H
#include "mat.h"
#endif

#include "Ptr.h"

/**
 * @brief Is matrix is a vector.
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
 * @brief   Matrix template.<br/>
 *          Core data structure
 *          
 * @author  Kaveh Vahedipour
 * @date    Mar 2010
 */
template <class T, paradigm P=SHM, const bool& b = TypeTraits<T>::Valid>
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
	inline
    Matrix              () {

	    //TypeTraits<T>::Validate();

        _dim.resize(1,1);
        _res.resize(1,1.0);
        
        Allocate();
        
    }
	
	
    /**
     * @brief           Construct matrix with dimension array
     *
     * @param  dim      All dimensions
     */
	inline explicit
    Matrix              (const std::vector<size_t>& dim) {
		
	    //TypeTraits<T>::Validate();

	    assert(!dim.empty() &&
	    		std::find(dim.begin(),dim.end(),size_t(0))==dim.end());

	    size_t n = 1, i = 0, ds = dim.size();

		_dim = dim;
		_res.resize(ds,1.0);

        Allocate();
        
	}
	
	
    
    /**
     * @brief           Construct matrix with aligned dimension vector
     *
     * @param  dim      All dimensions
     */
	inline
    Matrix              (const container<size_t>& dim) {
		
	    //TypeTraits<T>::Validate();

	    size_t ds = dim.size();
	    assert(ds &&
	    	   std::find(dim.begin(),dim.end(),size_t(0))==dim.end());

	    size_t n = 1, i = 0;
		
		_dim.resize(ds);
		std::copy (dim.begin(),dim.end(),_dim.begin());
		_res.resize(ds,1.0);

        Allocate();
		
	}
	
	
    /**
     * @brief           Construct matrix with dimension and resolution arrays
     *
     * @param  dim      All 16 Dimensions
     * @param  res      All 16 Resolutions
     */
	inline explicit
    Matrix              (const std::vector<size_t>& dim, const std::vector<float>& res) {
		
	    //TypeTraits<T>::Validate();

	    assert(!dim.empty() &&
	    	    std::find(dim.begin(),dim.end(),size_t(0))==dim.end() &&
	    	    dim.size() == res.size());

	    size_t n = 1, i = 0, ds = dim.size();
		
		_dim = dim;
		_res = res;

        Allocate();
		
	}

    
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
    inline explicit
    Matrix (const size_t n) {

	    //TypeTraits<T>::Validate();

	    assert (n);

		_dim.resize(2,n);
	    _res.resize(2,1.0);

        Allocate();
        
	}



    
    
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
    inline
    Matrix              (const size_t m, const size_t n) {

        //TypeTraits<T>::Validate();
    	assert (m && n);

    	_dim.resize(2); _dim[0] = m; _dim[1] = n;
		_res.resize(2,1.0);

        Allocate();

    }

    
    
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
    inline
    Matrix (const size_t m, const size_t n, const size_t k) {

	    //TypeTraits<T>::Validate();

	    assert (m && n && k);

    	_dim.resize(3); _dim[0] = m; _dim[1] = n; _dim[2] = k;
		_res.resize(3,1.0);

        Allocate();

    }
    

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
                         const size_t ave = 1) {

    	//TypeTraits<T>::Validate();

		assert (col && lin && cha && set && eco && phs && rep && seg &&
				par && slc && ida && idb && idc && idd && ide && ave );

	    size_t nd = 16, n = nd, sd, i;

	    _dim.resize(n);
	    _dim[ 0] = col; _dim[ 1] = lin; _dim[ 2] = cha; _dim[ 3] = set;
	    _dim[ 4] = eco; _dim[ 5] = phs; _dim[ 6] = rep; _dim[ 7] = seg;
	    _dim[ 8] = par; _dim[ 9] = slc; _dim[10] = ida; _dim[11] = idb;
	    _dim[12] = idc; _dim[13] = idd; _dim[14] = ide; _dim[15] = ave;

        // Remove trailing singleton dimensions
        for (i = 0; i < nd; i++)
            if (_dim[i] == 1)
                n--;
            else
                n = nd;
        
        // Resize skeleton
        _dim.resize(n);
		_res.resize(n,1.0);

        Allocate();

	}

    

    


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
    inline
    Matrix             (const Matrix<T,P> &M) {

		if (this != &M) {

		    //TypeTraits<T>::Validate();

			_dim = M.Dim();
			_dsz = M.Dsz();
			_res = M.Res();
	        _M   = M.Container();

		}

	}


    /**
     * @brief           Delete array containing data.
     */
	inline
	~Matrix() {};



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
    inline T
    operator[]          (const size_t p) const {
        assert(p <  Size());
        return _M[p];
    }
    
    
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
    inline T&
    operator[] (const size_t p) {
        assert(p <  Size());
        return _M[p];
    }


    
    /**
     * @brief           Get pointer to memory starting at p-th (default:0) element.
     *  
     * @param  p        Position
     *
     * @return          Data 
     */
    inline const T*            
    Memory             (const size_t p = 0)  const {
        assert (p < Size());
        return _M.memory(p);
    }

    
    /**
     * @brief           Data container (lhs)
     *  
     * @return          Data container
     */
    inline container<T>&
    Container           ()  {
        return _M;
    }

    
    /**
     * @brief           Data container (rhs)
     *  
     * @return          Data container
     */
    inline container<T>
    Container           ()  const {
        return _M;
    }


    /**
     * @brief           Container iterator to first element (lhs)
     *
     * @return          Container iterator
     */
    inline typename container<T>::iterator

    Begin               () {
    	return _M.begin ();
    }


    /**
     * @brief           Container const iterator to first element (rhs)
     *
     * @return          Container const iterator
     */
    inline typename container<T>::const_iterator
    Begin               ()  const {
    	return _M.begin ();
    }

    
    /**
     * @brief           Container iterator to last element (lhs)
     *
     * @return          Container iterator
     */
    inline typename container<T>::iterator
    End                 () {
    	return _M.end ();
    }


    /**
     * @brief           Container const iterator to last element (rhs)
     *
     * @return          Container const iterator
     */
    inline typename container<T>::const_iterator
    End                 ()  const {
    	return _M.end ();
    }


    /**
     * @brief           Element at position p (rhs)
     *  
     * @param  p        Position
     * @return          Value at _M[p]
     */
    inline T            
    At                  (const size_t p) const {
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
    At                  (const size_t pos) {
        assert (pos < Size());
        return _M[pos];
    }


    
    /**
     * @brief           Get element in (first) slice (rhs)
     *  
     * @param  x        Column
     * @param  y        Line
     *
     * @return          Value
     */
    inline T            
    At                  (const size_t x, const size_t y) const {

    	assert ((!x || x < _dim[0]) &&
    			(!y || y < _dim[1]));

        return _M[x + _dim[0]*y];

    }

    

    /**
     * @brief            Reference to element in (first) slice (rhs)
     *  
     * @param  x         Column
     * @param  y         Line
     *
     * @return           Reference
     */
    inline T&           
    At                  (const size_t x, const size_t y) {

    	assert ((!x || x < _dim[0]) && (!y || y < _dim[1]));

        return _M[x + _dim[0]*y];

    }

    
    /**
     * @brief          Get element in (first) volume (lhs)
     *  
     * @param  x       Column
     * @param  y       Line
     * @param  z       Slice
     *
     * @return         Value
     */
    inline T            
    At                   (const size_t x, const size_t y, const size_t z) const {

    	assert ((!x || x < _dim[0]) && (!y || y < _dim[1]) && (!z || z < _dim[2]));

        return _M[x + _dsz[1]*y + _dsz[2]*z];

    }
    
    

    /**
     * @brief            Reference to element in (first) volume (lhs)
     *  
     * @param  x         Column
     * @param  y         Line
     * @param  z         Slice
     *
     * @return           Reference
     */
    inline T&            
    At                   (const size_t x, const size_t y, const size_t z) {

    	assert ((!x || x < _dim[0]) && (!y || y < _dim[1]) && (!z || z < _dim[2]));

        return _M[x + _dsz[1]*y + _dsz[2]*z];

    }
    
    

    /**
     * @brief            Get value at position (lhs)
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

    	assert ((!col || col < _dim[ 0])&&
    	        (!lin || lin < _dim[ 1])&&
    	        (!cha || cha < _dim[ 2])&&
    	        (!set || set < _dim[ 3])&&
    	        (!eco || eco < _dim[ 4])&&
    	        (!phs || phs < _dim[ 5])&&
    	        (!rep || rep < _dim[ 6])&&
    	        (!seg || seg < _dim[ 7])&&
    	        (!par || par < _dim[ 8])&&
    	        (!slc || slc < _dim[ 9])&&
    	        (!ida || ida < _dim[10])&&
    	        (!idb || idb < _dim[11])&&
    	        (!idc || idc < _dim[12])&&
    	        (!idd || idd < _dim[13])&&
    	        (!ide || ide < _dim[14])&&
    	        (!ave || ave < _dim[15]));

        return _M [col + lin*_dsz[ 1] + cha*_dsz[ 2] + set*_dsz[ 3] + eco*_dsz[ 4] + phs*_dsz[ 5] + rep*_dsz[ 6] +
   	               seg*_dsz[ 7] + par*_dsz[ 8] + slc*_dsz[ 9] + ida*_dsz[10] + idb*_dsz[11] + idc*_dsz[12] +
   	               idd*_dsz[13] + ide*_dsz[14] + ave*_dsz[15]];

    }
    
    
    /**
     * @brief            Reference to value in volume (rhs)
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

    	assert ((!col || col < _dim[ 0])&&
    	       (!lin || lin < _dim[ 1])&&
    	       (!cha || cha < _dim[ 2])&&
    	       (!set || set < _dim[ 3])&&
    	       (!eco || eco < _dim[ 4])&&
    	       (!phs || phs < _dim[ 5])&&
    	       (!rep || rep < _dim[ 6])&&
    	       (!seg || seg < _dim[ 7])&&
    	       (!par || par < _dim[ 8])&&
    	       (!slc || slc < _dim[ 9])&&
    	       (!ida || ida < _dim[10])&&
    	       (!idb || idb < _dim[11])&&
    	       (!idc || idc < _dim[12])&&
    	       (!idd || idd < _dim[13])&&
    	       (!ide || ide < _dim[14])&&
    	       (!ave || ave < _dim[15]));

    	return _M [col + lin*_dsz[ 1] + cha*_dsz[ 2] + set*_dsz[ 3] + eco*_dsz[ 4] + phs*_dsz[ 5] + rep*_dsz[ 6] +
    	                 seg*_dsz[ 7] + par*_dsz[ 8] + slc*_dsz[ 9] + ida*_dsz[10] + idb*_dsz[11] + idc*_dsz[12] +
    	                 idd*_dsz[13] + ide*_dsz[14] + ave*_dsz[15]];

    }
    

    /**
     * @brief          Cast operator
     *
     * @return         Casted copy
     */
    template<class S>
    inline operator Matrix<S,P> () const {

		Matrix<S,P> m (_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			m[i] = (S)_M[i];

		return m;

	}




    /**
     * @brief           @see At(const size_t)
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    inline T
    operator()          (const size_t p) const {
        return this->At(p);
    }

    
    /**
     * @brief           Get value of pth element of repository.
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    inline T&
    operator()          (const size_t p) {
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
    operator()          (const size_t x, const size_t y) const {
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
    operator()           (const size_t x, const size_t y) {
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
    operator()           (const size_t x, const size_t y, const size_t z) const {
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
    operator()           (const size_t x, const size_t y, const size_t z) {
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
        return this->At (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave);
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
    inline Matrix<T,P>&
    operator=           (const Matrix<T,P>& M) {

        if (this != &M) {
        	_dim = M.Dim();
        	_dsz = M.Dsz();
        	_res = M.Res();
            _M  = M.Container();
        }

        return *this;

    }


    /**
     * @brief           Assignment data
     *
     * @param  v        Data vector (size must match numel(M)).
     */
    inline Matrix<T,P>&
    operator=           (const container<T>& v) {
        
    	assert (_M.size() == v.size());

        if (&_M != &v)
            _M = v;

        return *this;

    }
    
    
    
    /**
     * @brief           Assignment operator. Sets all elements s.
     *
     * @param  s        The assigned scalar.
     */
    inline Matrix<T,P>&
    operator=           (const T s) {

        T t = T(s);

#if defined USE_VALARRAY
        _M = t;        
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
            _M[i] = t;
#endif
        
        return *this;
    }
    
    
    /**
     * @brief           Matrix product. i.e. this * M.
     *
     * @param  M        The factor.
     */
    inline Matrix<T,P>
    operator->*         (const Matrix<T,P>& M) const {
        return this->prod(M);
    }

    /**
     * @brief           Elementwise substruction of two matrices
     *
     * @param  M        Matrix substruent.
     */
    inline Matrix<T,P>
    operator-           (const Matrix<T,P>& M) const {

		Matrix<T,P> res = *this;
		return res -= M;
        
    }

    
    /**
     * @brief           Elementwise substruction of two matrices
     *
     * @param  M        Matrix substruent.
     */
    template <class S>
    inline Matrix<T,P>
    operator-           (const Matrix<S,P>& M) const {

		Matrix<T,P> res = *this;
		return res -= M;

    }


    /**
     * @brief           Elementwise subtraction all elements by a scalar
     *
     * @param  s        Scalar substruent.
     */
    template <class S>
    inline Matrix<T,P>
    operator-           (const S s) const {

        T t = T(s);
        
        Matrix<T,P> res = *this;

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] -= t;
        
        return res;

    }

    
    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator-=         (const Matrix<T,P>& M) {

        size_t i;

        assert (Size() == M.Size());

#if defined USE_VALARRAY
        _M -= M.Container();
#elif defined EXPLICIT_SIMD
        if (fp_type(_M[0]))
        	SSE::binary<T>(_M, M.Container(), SSE::sub<T>(), _M);
        else
        	for (size_t i = 0; i < Size(); ++i)
        		_M[i] -= M[i];
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
            _M[i] -= M[i];
#endif

        return *this;

    }


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>&
    operator-=          (const Matrix<S,P>& M) {

        assert (Size() == M.Size());

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			_M[i] -= M[i];

        return *this;

    }

    /**
     * @brief           ELementwise substration with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
     * @return          Result
     */
    template <class S >
    inline Matrix<T,P>&
    operator-=          (const S s) {

		T t = T (s);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			_M[i] -= t;

		return *this;

    }

    /**
     * @brief           Unary minus (additive inverse)
     *
     * @return          Negation
     */
    inline Matrix<T,P>
    operator-           () const {

        Matrix<T,P> res = *this;

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] = -res[i];

        return res;

    }


    /**
     * @brief           Unary plus
     *
     * @return          Identity
     */
    inline Matrix<T,P>
    operator+           () const {

        return *this;

    }
    
    
    /**
     * @brief           Elementwise addition. i.e. this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>
    operator+          (const Matrix<T,P> &M) const {
        
		Matrix<T,P> res = *this;
		return res += M;
        
    }
    
    
    /**
     * @brief           Elementwise addition of two matrices
     *
     * @param  M        Matrix additive.
     */
    template <class S>
    inline const Matrix<T,P>
    operator+          (const Matrix<S,P>& M) const {
        
		Matrix<T,P> res = *this;
		return res += M;
		
    }
    
    
    /**
     * @brief           Elementwise addition iof all elements with a scalar
     *
     * @param  s        Scalar additive.
     */
    template <class S>
    inline Matrix<T,P>
    operator+           (const S s) const {

        Matrix<T,P> res = *this;
    	T t = T(s);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] += t;

        return res;

    }

    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator+=         (const Matrix<T,P>& M) {

        assert (Size() == M.Size());

#if defined USE_VALARRAY
        _M += M.Container();        
#elif defined EXPLICIT_SIMD
        if (fp_type(_M[0]))
        	SSE::binary<T>(_M, M.Container(), SSE::add<T>(), _M);
        else
        	for (size_t i = 0; i < Size(); ++i)
        		_M[i] += M[i];
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
            _M[i] += M[i];
#endif
        
        return *this;
        
    }

    
/**
 * @brief           ELementwise multiplication and assignment operator. i.e. this = m.
 *
 * @param  M        Added matrix.
 * @return          Result
 */
    template <class S>
    inline Matrix<T,P>&
    operator+=          (const Matrix<S,P>& M) {

        assert (Size() == M.Size());

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			_M[i] += M[i];

    	return *this;

    }


    /**
     * @brief           ELementwise addition with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
     * @return          Result
     */
    template <class S >
    inline Matrix<T,P>&
    operator+=          (const S s) {

    	T t = T (s);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			_M[i] += t;

        return *this;

    }

    
    
    /**
     * @brief           Transposition / Complex conjugation. i.e. this'.
     *
     * @return          Matrix::tr()
     */
    inline Matrix<T,P>
    operator!           () const {

        Matrix<T,P> res (_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < _dim[1]; ++i)
			for (size_t j = 0; j < _dim[0]; j++)
				res(i,j) = this->At(j,i);

        return res;

    }

    
    
    /**
     * @brief           Return a matrix with result[i] = (m[i] ? this[i] : 0).
     *
     * @param  M        The operand
     * @return          Cross-section or zero
     */
    inline Matrix<T,P>
    operator&           (const Matrix<cbool>& M) const ;
    
    
     /**
     * @brief           Scalar inequality. result[i] = (this[i] != m). i.e. this ~= m
     *
     * @param  s        Comparing scalar.
     * @return          Matrix of false where elements are equal s and true else.
     */
    inline Matrix<cbool>
    operator!=          (const T s) const {

        Matrix<cbool> res(_dim);
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = (_M[i] != s) ? 1 : 0;
        return res;

    }

    
    
    /**
     * @brief           Scalar greater comaprison, result[i] = (this[i] > m). i.e. this > m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator>           (const T s) const {

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::greater(_M[i], s);

        return res;

    }

    
    
    /**
     * @brief           Scalar greater or equal comparison. result[i] = (this[i] >= m). i.e. this >= m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator>=          (const T s) const {

		Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::greater_or_equal(_M[i], s);

        return res;

	}

    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] <= m). i.e. this <= m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator<=          (const T s) const {

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::less_or_equal(_M[i], s);

        return res;

    }

    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] < m). i.e. this < m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator<           (const T s) const {

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::less(_M[i], s);

        return res;

    }


    /**
     * @brief           Elementwise equality, result[i] = (this[i] == m[i]). i.e. this == m
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator==          (const Matrix<T,P>& M) const {

        assert (Size() == M.Size());

        Matrix<cbool> res(_dim);
#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] = (_M[i] == M[i]) ? 1 : 0;

        return res;

    }


    /**
	 * @brief           Elementwise equality, result[i] = (this[i] == m[i]). i.e. this == m
	 *
	 * @param  M        Comparing matrix.
	 * @return          Hit list
	 */
    template<class S>
	inline Matrix<cbool>
	operator==          (const Matrix<S,P>& M) const {

        assert (Size() == M.Size());

		Matrix<cbool> res (_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] = (_M[i] == (T)M[i]) ? 1 : 0;

		return res;

     }

    /**
     * @brief           Scalar equality. result[i] = (this[i] == m).
     *
     * @param  s        Comparing scalar.
     * @return          Matrix of true where elements are equal s and false else.
     */
    inline Matrix<cbool>
    operator==          (const T s) const {

    	T t = (T) s;

        Matrix<cbool> res (_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] =  (_M[i] == s) ? 1 : 0;

        return res;

    }


    
    /**
     * @brief           Elementwise equality, result[i] = (this[i] != m[i]). i.e. this ~= m
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator!=          (const Matrix<T,P>& M) const {

        assert (Size() == M.Size());

        Matrix<cbool> res(_dim,_res);
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = (_M[i]!= M[i]) ? 1 : 0;
        return res;

    }

    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] >= m[i]). i.e. this >= m
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator>=          (const Matrix<T,P>& M) const {

        assert (Size() == M.Size());

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::greater_or_equal(_M[i], M[i]);

        return res;

    }

    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] <= m[i]). i.e. this <= m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator<=          (const Matrix<T,P>& M) const {

        assert (Size() == M.Size());

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::less_or_equal(_M[i], M[i]);

        return res;

    }

    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] > m[i]). i.e. this > m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator>           (const Matrix<T,P>& M) const {

        assert (Size() == M.Size());

        Matrix<cbool> res(_dim,_res);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::greater(_M[i], M[i]);

        return res;

    }

    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] < m[i]). i.e. this < m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator<           (const Matrix<T,P>& M) const {

        assert (Size() == M.Size());

        Matrix<cbool> res(_dim,_res);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::less(_M[i], M[i]);

        return res;

    }

    
    /**
     * @brief           Matrix comparison, result[i] = (m[i] || this[i] ? 1 : 0). i.e. this | m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator||          (const Matrix<T,P>& M) const {

        assert (Size() == M.Size());

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::logical_or(_M[i], M[i]);

        return res;

    }

    
    
    /**
     * @brief           Matrix comparison, result[i] = (m[i] && this[i] ? 1 : 0). i.e. this & m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator&&          (const Matrix<T,P>& M) const {

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::logical_and(_M[i], M[i]);

        return res;

    }


    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Power.
     * @return          Result
     */
    inline Matrix<T,P>
    operator^           (const float p) const {

    	Matrix<T,P> res = *this;

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] = (p == 0) ? T(1) : pow(res[i],  p);

        return res;

    }


    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Power.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator^=          (const float p) {

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			_M[i] = pow(_M[i],  p);

        return *this;

    }
    

	#include "Operators.hpp"
    
    

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

#ifdef HAVE_SCALAPACK

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
     * @brief           Get resolution of a given dimension.
     *
     * @param   i       Dimension
     * @return          Resolution .
     */
    inline float
    Res                 (const size_t i) const {
        assert (i < _res.size());
        return _res[i];
    }


    /**
     * @brief           Resolution of a given dimension.
     *
     * @param   i       Dimension
     * @return          Resolution
     */
    inline float&
    Res                 (const size_t i)       {
        assert (i < _res.size());
        return _res[i];
    }



    /**
     * @brief           Resolution array
     *
     * @return          All resolutions
     */
    inline std::vector<float>
    Res                 () const {
        return _res;
    }



    /**
     * @brief           Get size a given dimension.
     *
     * @param   i       Dimension
     * @return          Dimension
     */
    inline size_t
    Dim                 (const size_t i)  const {
        return (i < _dim.size()) ? _dim[i]: 1;
    }


    /**
     * @brief           Get dimension array
     *
     * @return          All dimensions
     */
    inline std::vector<size_t>
    Dim                 ()                  const {
        return _dim;
    }


    inline size_t
    NDim() const {
    	return _dim.size();
    }


    /**
     * @brief           Get dimension sizes
     *
     * @return          All dimensions
     */
    inline std::vector<size_t>
    Dsz                 ()                  const {
        return _dsz;
    }


    /**
     * @brief           Purge data and free RAM.
     */
    inline void
    Clear               ()                                      {

    	_dim.resize(1,1);
        _dsz.resize(1,1);
        _res.resize(1,1.0);
        _M.resize(1);

    }


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
    inline const char*
    GetClassName        () const { 
        return _name.c_str(); 
    }


    /**
     * @brief           Who are we?
     *
     * @return          Class name
     */ 
    inline void
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
    inline Matrix<T,P>
    prod                (const Matrix<T,P> &M, const char transa = 'N', const char transb = 'N') const {
        return gemm (*this, M, transa, transb);
    }


    /**
     * @brief           Complex conjugate left and multiply with right.
     *
     * @param   M       Factor
     * @return          Product of conj(this) and M.
     */
    inline Matrix<T,P>
    prodt               (const Matrix<T,P> &M) const {
        return gemm (*this, M, 'C');
    }


    /**
     * @brief           Scalar product (complex: conjugate first vector) using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTC and XDOT
     *
     * @param  M        Factor
     * @return          Scalar product
     */
    inline T
    dotc (const Matrix<T,P>& M) const  {
        return DOTC (*this, M);
    }


    /**
     * @brief           Scalar product using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTU and XDOT
     *
     * @param  M        Factor
     * @return          Scalar product
     */
    inline T
    dot (const Matrix<T,P>& M) const {

        return DOT  (*this, M);

    }

    
    
    /**
     * @brief           Number of elements
     *
     * @return          Size
     */
    inline size_t
    Size () const {        
        return _M.size();
    }

    //@}


protected:
	
    /**
     * @brief           Number of elements
     *
     * @return          Size
     */
    inline size_t
    DimProd () const {
        
        long size = 1;
        
        for (size_t i = 0; i < _dim.size(); ++i)
            size *= _dim[i];
        
        return size;
        
    }

    inline void
    Allocate () {

        size_t ds = _dim.size(), i;
		_dsz.resize(ds,1);
	    for (i = 1; i < ds; ++i)
	        _dsz[i] = _dsz[i-1]*_dim[i-1];
        
        _M = container<T>(DimProd());
        
    }

    // Structure
    std::vector<size_t> _dim; /// Dimensions
    std::vector<size_t> _dsz; /// Dimension size.
    std::vector<float>  _res; /// Resolutions
    
    //Data
    container<T>    _M; /// Data container

    // Name
    std::string         _name; /// Name
    
#ifdef HAVE_SCALAPACK
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

};


#ifdef PARC_MODULE_NAME

template <class T, paradigm P> long 
Matrix<T,P>::Import     (const IceAs* ias, const size_t pos) {
    
    ICE_SET_FN("Matrix<T,P>::Import(IceAs, long)")
        
    int  i    = 0;
    long size = 1;
    
    for (i = 0; i < MAX_ICE_DIM; ++i)
        size *= (ias->getLen(IceDim(i)) <= 1) ? 1 : ias->getLen(IceDim(i));
    
    T* data = (T*) ias->calcSplObjStartAddr() ;
    
    for (i = 0; i < size; ++i, ++data)
        _M[i+pos] = *data;
    
    return size;
    
}


template <class T, paradigm P> long 
Matrix<T,P>::Import(const IceAs* ias) {
    
    ICE_SET_FN("Matrix<T,P>::Import(IceAs)")
        
    int i;
    
    for (i = 0; i < MAX_ICE_DIM; ++i)
        _dim[i] = (ias->getLen(IceDim(i)) <= 1) ? 1 : ias->getLen(IceDim(i));
    
    _M = new T[Size()]();
    nb_alloc = 1;
    
    T* data = (T*) ias->calcSplObjStartAddr() ;
    
    for (i = 0; i < Size(); ++i, data++)
        _M[i] = *data;
    
    return Size();
    
}


template <class T, paradigm P> long 
Matrix<T,P>::Export (IceAs* ias) const {
    
    ICE_SET_FN("Matrix<T,P>::Export(IceAs)")
        
    T* data = (T*) ias->calcSplObjStartAddr() ;
    
    for (int i = 0; i < Size(); ++i, data++)
        *data = _M[i];
    
    return Size();
    
}


template <class T, paradigm P> long
Matrix<T,P>::Export (IceAs* ias, const size_t pos) const {

    ICE_SET_FN("Matrix<T,P>::Export(IceAs, long)")
        
        int  i    = 0;
    long size = 1;
    
    for (i = 0; i < MSX_ICE_DIM; ++i) {
        size *= (ias->getLen(IceDim(i)) <= 1) ? 1 : ias->getLen(IceDim(i));
    }
    
    T* data = (T*) ias->calcSplObjStartAddr() ;
    
    for (i = 0; i < size; ++i, data++)
        *data = _M[i+pos];
    
    return size;
    
}

#endif // ICE

#ifdef HAVE_SCALAPACK
#include "Grid.hpp"

template<> inline
Matrix<float,MPI>::Matrix (const size_t cols, const size_t rows) {

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
    //this->_M = VECTOR_CONSTR(float,Size());
	
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


