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
#include "Vector.hpp"

#ifdef EXPLICIT_SIMD
#    include "SIMD.hpp"
#endif

#include <iostream>
#include <memory>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <numeric>

#include <ostream>
#include <string>
#include <cstring>
#include <algorithm>
#include <utility>

#include <Assert.hpp>


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
template <class T, paradigm P=SHM
#if !defined(_MSC_VER) || _MSC_VER>1200
    ,const bool& b = TypeTraits<T>::Supported
#endif
>
class Matrix {
	
	
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

        _dim.resize(1,1);
        _res.resize(1,1.0);
        
        Allocate();
        
    }
	
	
    /**
     * @brief           Construct matrix with aligned dimension vector
     *
     * @param  dim      All dimensions
     */
	inline
    Matrix              (const Vector<size_t>& dim) {

	    size_t ds = dim.size();
	    assert(ds &&
	    	   std::find(dim.begin(),dim.end(),size_t(0))==dim.end());

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
    Matrix              (const Vector<size_t>& dim, const Vector<float>& res) {
		
	    assert(!dim.Empty() &&
	    	    std::find(dim.begin(),dim.end(),size_t(0))==dim.end() &&
	    	    dim.size() == res.size());

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
    Matrix              (const size_t col,     const size_t lin,     const size_t cha,     const size_t set,
                         const size_t eco = 1, const size_t phs = 1, const size_t rep = 1, const size_t seg = 1,
                         const size_t par = 1, const size_t slc = 1, const size_t ida = 1, const size_t idb = 1,
                         const size_t idc = 1, const size_t idd = 1, const size_t ide = 1, const size_t ave = 1) {

		assert (col && lin && cha && set && eco && phs && rep && seg &&
				par && slc && ida && idb && idc && idd && ide && ave );

	    size_t nd = 16, n = nd, i;

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


#ifdef HAVE_CXX11_RVALUE_REFERENCES
    /**
     * @brief           Default copy constructor
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
    	if (this != &M)
    		*this = M;
    }
    /**
     * @brief           Default move constructor
     *
     * Usage:
     * @code{.cpp}
     *   Matrix<cxfl> m (n); // Copy n into m
     * @endcode
     *
     * @param  M        Right hand side
     */
    inline
    Matrix             (Matrix<T,P>&& M) {
    	if (this != &M)
    		*this = M;
    }
#endif

    //@}



    /**
     * @name            Import export functions for ICE access specifiers.<br/>
     *                  Ice access specifiers can be handled in one of the following ways.<br/>
     *                  It is crucial to understand that the
     */

    //{@

    // Only if compiled within IDEA we know of access specifiers.
#ifdef PARC_MODULE_NAME
    #include "ICE.hpp"
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
     * @brief           Get reference to p-th element.
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


	//#include "AssignmentHandler.hpp"
    
    //inline AssignmentHandler<T> operator[] (const int, const int) {

    //}

    /**
     * @brief           Get pointer to memory starting at p-th (default:0) element.
     *  
     * @param  p        Position
     *
     * @return          Data 
     */
    inline const T*            
    Ptr             (const size_t p = 0)  const {
        assert (p < Size());
        return _M.ptr(p);
    }

    /**
     * @brief           Get pointer to memory starting at p-th (default:0) element.
     *
     * @param  p        Position
     *
     * @return          Data
     */
    inline T*
    Ptr             (const size_t p = 0) {
        assert (p < Size());
        return _M.ptr(p);
    }


    
    /**
     * @brief           Data container (lhs)
     *  
     * @return          Data container
     */
    inline Vector<T>&
    Container           ()  {
        return _M;
    }

    
    /**
     * @brief           Data container (rhs)
     *  
     * @return          Data container
     */
    inline Vector<T>
    Container           ()  const {
        return _M;
    }


    /**
     * @brief           Container iterator to first element (lhs)
     *
     * @return          Container iterator
     */
    inline typename Vector<T>::iterator

    Begin               () {
    	return _M.begin ();
    }


    /**
     * @brief           Container const iterator to first element (rhs)
     *
     * @return          Container const iterator
     */
    inline typename Vector<T>::const_iterator
    Begin               ()  const {
    	return _M.begin ();
    }

    
    /**
     * @brief           Container iterator to last element (lhs)
     *
     * @return          Container iterator
     */
    inline typename Vector<T>::iterator
    End                 () {
    	return _M.end ();
    }


    /**
     * @brief           Container const iterator to last element (rhs)
     *
     * @return          Container const iterator
     */
    inline typename Vector<T>::const_iterator
    End                 ()  const {
    	return _M.end ();
    }


    /**
     * @brief           Element at position p (rhs)
     *  
     * @param  p        Position
     * @return          Value at _M[p]
     */
    inline const T&
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
    inline const T& At (const size_t x, const size_t y) const {
    	assert ((!x || x < _dim[0]) && (!y || y < _dim[1]));
        return _M[x + _dim[0]*y];
    }
    inline T&       At (const size_t x, const size_t y) {
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
    inline const T& At (const size_t x, const size_t y, const size_t z) const {
    	assert ((!x || x < _dim[0]) && (!y || y < _dim[1]) && (!z || z < _dim[2]));
        return _M[x + _dsz[1]*y + _dsz[2]*z];
    }
    inline T&       At (const size_t x, const size_t y, const size_t z) {
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
    inline const T&
    At                   (const size_t x, const size_t y, const size_t z, const size_t w) const {
    	assert ((!x || x < _dim[0]) && (!y || y < _dim[1]) && (!z || z < _dim[2]) && (!w || w < _dim[3]));
        return _M[x + _dsz[1]*y + _dsz[2]*z + _dsz[3]*w];
    }
    inline T&
    At                   (const size_t x, const size_t y, const size_t z, const size_t w) {
    	assert ((!x || x < _dim[0]) && (!y || y < _dim[1]) && (!z || z < _dim[2]) && (!w || w < _dim[3]));
        return _M[x + _dsz[1]*y + _dsz[2]*z + _dsz[3]*w];
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
    inline const T&
    At                   (const size_t col,     const size_t lin,     const size_t cha,     const size_t set,
                          const size_t eco,     const size_t phs = 0, const size_t rep = 0, const size_t seg = 0,
                          const size_t par = 0, const size_t slc = 0, const size_t ida = 0, const size_t idb = 0,
                          const size_t idc = 0, const size_t idd = 0, const size_t ide = 0, const size_t ave = 0) const {
    	assert ((!col || col < _dim[ 0])&& (!lin || lin < _dim[ 1])&& (!cha || cha < _dim[ 2])&& (!set || set < _dim[ 3])&&
    	        (!eco || eco < _dim[ 4])&& (!phs || phs < _dim[ 5])&& (!rep || rep < _dim[ 6])&& (!seg || seg < _dim[ 7])&&
    	        (!par || par < _dim[ 8])&& (!slc || slc < _dim[ 9])&& (!ida || ida < _dim[10])&& (!idb || idb < _dim[11])&&
    	        (!idc || idc < _dim[12])&& (!idd || idd < _dim[13])&& (!ide || ide < _dim[14])&& (!ave || ave < _dim[15]));
        return _M [col + lin*_dsz[ 1] + cha*_dsz[ 2] + set*_dsz[ 3] + eco*_dsz[ 4] + phs*_dsz[ 5] + rep*_dsz[ 6] +
   	               seg*_dsz[ 7] + par*_dsz[ 8] + slc*_dsz[ 9] + ida*_dsz[10] + idb*_dsz[11] + idc*_dsz[12] +
   	               idd*_dsz[13] + ide*_dsz[14] + ave*_dsz[15]];
    }
    inline T&
    At                   (const size_t col,     const size_t lin,     const size_t cha,     const size_t set,
                          const size_t eco,     const size_t phs = 0, const size_t rep = 0, const size_t seg = 0,
                          const size_t par = 0, const size_t slc = 0, const size_t ida = 0, const size_t idb = 0,
                          const size_t idc = 0, const size_t idd = 0, const size_t ide = 0, const size_t ave = 0) {
    	assert ((!col || col < _dim[ 0])&& (!lin || lin < _dim[ 1])&& (!cha || cha < _dim[ 2])&& (!set || set < _dim[ 3])&&
    	        (!eco || eco < _dim[ 4])&& (!phs || phs < _dim[ 5])&& (!rep || rep < _dim[ 6])&& (!seg || seg < _dim[ 7])&&
    	        (!par || par < _dim[ 8])&& (!slc || slc < _dim[ 9])&& (!ida || ida < _dim[10])&& (!idb || idb < _dim[11])&&
    	        (!idc || idc < _dim[12])&& (!idd || idd < _dim[13])&& (!ide || ide < _dim[14])&& (!ave || ave < _dim[15]));
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
		for (size_t i = 0; i < Size(); ++i)
			m[i] = (S)_M[i];
		return m;
	}

    /**
     * @brief Deliver range of values with indices
     *
     * @param indices List of indices
     * @return        Matrix containing values at indices
     */
    inline Matrix<T> operator() (const Vector<size_t>& indices, size_t col = 0) const {
    	size_t indices_size = indices.size();
    	assert (indices_size);
    	assert (indices_size < Size());
    	assert (mmax(indices) < Size());
    	Matrix<T> ret (indices_size,1);
    	for (size_t i = 0; i < indices_size; ++i)
    		ret[i] = At(indices[i],col);
    	return ret;
    }

    /**
     * @brief Deliver range of values with indices
     *
     * @param indices List of indices
     * @return        Matrix containing values at indices
     */
    inline Matrix<T> operator() (size_t row, const Vector<size_t>& indices) const {
    	size_t indices_size = indices.size();
    	assert (indices_size);
    	assert (indices_size < Size());
    	assert (mmax(indices) < Size());
    	Matrix<T> ret (1,indices_size);
    	for (size_t i = 0; i < indices_size; ++i)
    		ret[i] = At(row,indices[i]);
    	return ret;
    }

    /**
     * @brief Deliver range of values with indices
     *
     * @param indices List of indices
     * @return        Matrix containing values at indices
     */
    inline Matrix<T> operator() (const Vector<size_t>& row_inds, const Vector<size_t>& col_inds) const {
    	size_t row_inds_size = row_inds.size();
    	size_t col_inds_size = col_inds.size();
    	assert (row_inds_size);
    	assert (col_inds_size);
    	assert (row_inds_size < _dim[0]);
    	assert (col_inds_size < _dim[1]);
    	assert (mmax(row_inds) < _dim[0]);
    	assert (mmax(col_inds) < _dim[1]);
    	Matrix<T> ret (row_inds_size, col_inds_size);
    	for (size_t j = 0; j < col_inds_size; ++j)
    		for (size_t i = 0; i < row_inds_size; ++i)
    			ret(i,j) = At(row_inds[i],col_inds[i]);
    	return ret;
    }




    /**
     * @brief           @see At(const size_t)
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    inline const T&
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
    inline const T&
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
    inline const T&
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
     * @brief            Get value in volume
     *
     * @param  x         Column
     * @param  y         Line
     * @param  z         Slice
     *
     * @return           Value
     */
    inline const T&
    operator()           (const size_t x, const size_t y, const size_t z, const size_t w) const {
        return this->At(x,y,z,w);
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
    operator()           (const size_t x, const size_t y, const size_t z, const size_t w) {
        return this->At(x,y,z,w);
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
    operator()     (const size_t col,     const size_t lin,     const size_t cha,     const size_t set,
					const size_t eco,     const size_t phs = 0, const size_t rep = 0, const size_t seg = 0,
					const size_t par = 0, const size_t slc = 0, const size_t ida = 0, const size_t idb = 0,
					const size_t idc = 0, const size_t idd = 0, const size_t ide = 0, const size_t ave = 0) {

		assert ((!col || col < _dim[ 0])&& (!lin || lin < _dim[ 1])&& (!cha || cha < _dim[ 2])&& (!set || set < _dim[ 3])&&
		  (!eco || eco < _dim[ 4])&& (!phs || phs < _dim[ 5])&& (!rep || rep < _dim[ 6])&& (!seg || seg < _dim[ 7])&&
		  (!par || par < _dim[ 8])&& (!slc || slc < _dim[ 9])&& (!ida || ida < _dim[10])&& (!idb || idb < _dim[11])&&
		  (!idc || idc < _dim[12])&& (!idd || idd < _dim[13])&& (!ide || ide < _dim[14])&& (!ave || ave < _dim[15]));

		return _M [col + lin*_dsz[ 1] + cha*_dsz[ 2] + set*_dsz[ 3] + eco*_dsz[ 4] + phs*_dsz[ 5] + rep*_dsz[ 6] +
				seg*_dsz[ 7] + par*_dsz[ 8] + slc*_dsz[ 9] + ida*_dsz[10] + idb*_dsz[11] + idc*_dsz[12] +
				idd*_dsz[13] + ide*_dsz[14] + ave*_dsz[15]];

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
    inline const T&
    operator()           (const size_t col,     const size_t lin,     const size_t cha,     const size_t set,
                          const size_t eco,     const size_t phs = 0, const size_t rep = 0, const size_t seg = 0,
                          const size_t par = 0, const size_t slc = 0, const size_t ida = 0, const size_t idb = 0,
                          const size_t idc = 0, const size_t idd = 0, const size_t ide = 0, const size_t ave = 0) const {

    	assert ((!col || col < _dim[ 0])&& (!lin || lin < _dim[ 1])&& (!cha || cha < _dim[ 2])&& (!set || set < _dim[ 3])&&
    	        (!eco || eco < _dim[ 4])&& (!phs || phs < _dim[ 5])&& (!rep || rep < _dim[ 6])&& (!seg || seg < _dim[ 7])&&
    	        (!par || par < _dim[ 8])&& (!slc || slc < _dim[ 9])&& (!ida || ida < _dim[10])&& (!idb || idb < _dim[11])&&
    	        (!idc || idc < _dim[12])&& (!idd || idd < _dim[13])&& (!ide || ide < _dim[14])&& (!ave || ave < _dim[15]));

        return _M [col + lin*_dsz[ 1] + cha*_dsz[ 2] + set*_dsz[ 3] + eco*_dsz[ 4] + phs*_dsz[ 5] + rep*_dsz[ 6] +
   	               seg*_dsz[ 7] + par*_dsz[ 8] + slc*_dsz[ 9] + ida*_dsz[10] + idb*_dsz[11] + idc*_dsz[12] +
   	               idd*_dsz[13] + ide*_dsz[14] + ave*_dsz[15]];

    }


    //@}
    

    
    
#ifndef NO_LAPACK
    /**
     * @brief           Matrix product. i.e. this * M.
     *
     * @param  M        The factor.
     */
    inline Matrix<T,P>
    operator->*         (const Matrix<T,P>& M) const;

    /**
     * @brief           Matrix Product.
     *
     * @param   M       The factor
     * @param   transa  Transpose ('T') / Conjugate transpose ('C') the left matrix. Default: No transposition'N'
     * @param   transb  Transpose ('T') / Conjugate transpose ('C') the right matrix. Default: No transposition 'N'
     * @return          Product of this and M.
     */
    Matrix<T,P>
    prod                (const Matrix<T,P> &M, const char transa = 'N', const char transb = 'N') const;



    /**
     * @brief           Complex conjugate left and multiply with right.
     *
     * @param   M       Factor
     * @return          Product of conj(this) and M.
     */
    Matrix<T,P>
    prodt               (const Matrix<T,P> &M) const;


    /**
     * @brief           Scalar product (complex: conjugate first vector) using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTC and XDOT
     *
     * @param  M        Factor
     * @return          Scalar product
     */
    inline T
    dotc (const Matrix<T,P>& M) const;


    /**
     * @brief           Scalar product using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTU and XDOT
     *
     * @param  M        Factor
     * @return          Scalar product
     */
    inline T
    dot (const Matrix<T,P>& M) const;
    
#endif
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
        return (_dim.size()>1) ? _dim[1] : 1;
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
    inline const Vector<float>&
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
    inline const Vector<size_t>&
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
    inline const Vector<size_t>&
    Dsz                 ()                  const {
        return _dsz;
    }


    /**
     * @brief           Purge data and free RAM.
     */
    inline void
    Clear               ()                                      {
    	_dim.Clear();
        _dsz.Clear();
        _res.Clear();
        _M.Clear();
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
     * @brief           Number of elements
     *
     * @return          Size
     */
    inline size_t
    Size () const {        
        return _M.size();
    }

    //@}


#include <functional>

/**
     * @name            Some operators
     *                  Operator definitions. Needs big expansion still.
     */

    //@{

#ifdef HAVE_CXX11_RVALUE_REFERENCES
	/**
	 * @brief  Move assignment operator
	 * @param  rhs  Right hand side reference
	 * @return      Left hand side
	 */
    inline Matrix<T,P>&
    operator= (Matrix<T,P>&& rhs) {
        _M    = std::move(rhs._M);
        _name = std::move(rhs._name);
        _res  = std::move(rhs._res);
        _dsz  = std::move(rhs._dsz);
        _dim  = std::move(rhs._dim);
        return *this;
    }
#endif
	/**
	 * @brief  Copy assignment operator
	 * @param  rhs  Right hand side
	 * @return      Left hand side
	 */
    inline Matrix<T,P>&
    operator= (const Matrix<T,P>& rhs) {
        _M    = rhs._M;
        _name = rhs._name;
        _res  = rhs._res;
        _dsz  = rhs._dsz;
        _dim  = rhs._dim;
        return *this;
    }

    /**
     * @brief           Assignment data
     *
     * @param  v        Data vector (size must match numel(M)).
     */
    inline Matrix<T,P>&
    operator=           (const Vector<T>& v) {

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
    inline const Matrix<T,P>&
    operator=           (const T s) {
        T t = T(s);
        std::fill(_M.begin(), _M.end(), t);
        return *this;
    }


    /**
     * @brief           Unary minus (additive inverse)
     *
     * @return          Negation
     */
    inline Matrix<T,P>
    operator-           () const {
        Matrix<T,P> res (_dim);
        std::transform (_M.begin(), _M.end(), res.Begin(), std::negate<T>());
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
     * @brief           Transposition / Complex conjugation. i.e. this'.
     *
     * @return          Matrix::tr()
     */
    inline Matrix<T,P>
    operator!           () const {
    	assert(_dim.size() ==2);
    	Matrix<T,P> res = *this;
		for (size_t i = 0; i < _dim[1]; ++i)
			for (size_t j = 0; j < i; j++)
				swapd(res(j,i),res(i,j));
        return res;
    }



    /**
     * @brief           Return a matrix with result[i] = (m[i] ? this[i] : 0).
     *
     * @param  M        The operand
     * @return          Cross-section or zero
     */
    //inline Matrix<T,P>
    //operator&           (const Matrix<cbool>& M) const ;


     /**
     * @brief           Scalar inequality. result[i] = (this[i] != m). i.e. this ~= m
     *
     * @param  s        Comparing scalar.
     * @return          Matrix of false where elements are equal s and true else.
     */
    inline Matrix<cbool>
    operator!=          (const T s) const {
        Matrix<cbool> res(_dim);
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
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::less(_M[i], s);
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
        op_assert (_dim == M.Dim(), *this, M);
        Matrix<cbool> res(_dim,_res);
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
    	op_assert (_dim == M.Dim(), *this, M);
        Matrix<cbool> res(_dim);
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
        op_assert (_dim == M.Dim(), *this, M);
        Matrix<cbool> res(_dim);
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
        op_assert (_dim == M.Dim(), *this, M);
        Matrix<cbool> res(_dim,_res);
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
        op_assert (_dim == M.Dim(), *this, M);
        Matrix<cbool> res(_dim,_res);
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
        Matrix<cbool> res(_dim);
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::logical_or(_M[i], M[i]);
        return res;
    }


    /**
     * @brief           Matrix comparison, result[i] = (m[i] || this[i] ? 1 : 0). i.e. this | m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<short>
    operator|          (const Matrix<T,P>& rhs) const {
        Matrix<cbool> ret(_dim);
        for (size_t i = 0; i < Size(); ++i)
        	ret[i] = (short)CompTraits<T>::logical_or(_M[i], ret[i]);
        return ret;
    }


    /**
     * @brief           Matrix comparison, result[i] = (m[i] && this[i] ? 1 : 0). i.e. this & m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator&&          (const Matrix<T,P>& M) const {
        op_assert (_dim == M.Dim(), *this, M);
        Matrix<cbool> res(_dim);
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
		for (size_t i = 0; i < Size(); ++i)
			res[i] = (p == 0) ? T(1) : TypeTraits<T>::Pow(res[i],  p);
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
		for (size_t i = 0; i < Size(); ++i)
			_M[i] = TypeTraits<T>::Pow(_M[i],  p);
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
        return res += s;
    }


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator+=         (const Matrix<T,P>& M) {
    	op_assert (_dim == M.Dim(), *this, M);
   		std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::plus<T>());
        return *this;
    }


/**
 * @brief           ELementwise multiplication and assignment operator. i.e. this = m.
 *
 * @param  M        Added matrix.
 * @return          Result
 */
    template <class S> inline Matrix<T,P>&
    operator+=         (const Matrix<S,P>& M) {
        op_assert (_dim == M.Dim(), *this, M);
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::plus<T>());
		return *this;
    }


    /**
     * @brief           ELementwise addition with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
     * @return          Result
     */
    template <class S > inline Matrix<T,P>&
    operator+=          (const S s) {
        std::transform (_M.begin(), _M.end(), _M.begin(), std::bind2nd(std::plus<T>(),(T)s));
        return *this;
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
    template <class S> inline Matrix<T,P>
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
		Matrix<T,P> res = *this;
		return res -= s;
    }



    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator-=         (const Matrix<T,P>& M) {
        op_assert (_dim == M.Dim(), *this, M);
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::minus<T>());
        return *this;
    }


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
     * @return          Result
     */
    template <class S> inline Matrix<T,P>&
    operator-=          (const Matrix<S,P>& M) {
        op_assert (_dim == M.Dim(), *this, M);
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::minus<T>());
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
        std::transform (_M.begin(), _M.end(), _M.begin(), std::bind2nd(std::minus<T>(),(T)s));
		return *this;
    }

    /**
     * @brief           Elementwise multiplication. i.e. this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>
    operator*          (const Matrix<T,P> &M) const {
		Matrix<T,P> res = *this;
		return res *= M;
    }


    /**
     * @brief           Elementwise multiplication. i.e. this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>
    operator*          (const Matrix<S,P> &M) const {
		Matrix<T,P> res = *this;
		return res *= M;
    }


    /**
     * @brief           Elementwise multiplication with a scalar. i.e. this * m.
     *
     * @param  s        Factor scalar
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>
    operator*          (const S s) const  {
        Matrix<T,P> res = *this;
        return res *= s;
    }



    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator*=         (const Matrix<T,P>& M) {
        op_assert (_dim == M.Dim(), *this, M);
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::multiplies<T>());
        return *this;
    }


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    template <class S> inline Matrix<T,P>&
    operator*=         (const Matrix<S,P>& M) {
        op_assert (_dim == M.Dim(), *this, M);
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::multiplies<T>());
        return *this;
    }


    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this *= s.
     *
     * @param  s        Factor scalar.
     * @return          Result
     */
    template <class S> inline Matrix<T,P>&
    operator*=         (const S s) {
        std::transform (_M.begin(), _M.end(), _M.begin(), std::bind2nd(std::multiplies<T>(),(T)s));
		return *this;
    }


    /**
     * @brief           Elementwise substruction of two matrices
     *
     * @param  M        Matrix substruent.
     */
    inline Matrix<T,P>
    operator/           (const Matrix<T,P>& M) const {
		Matrix<T,P> res = *this;
		return res /= M;
    }


    /**
     * @brief           Elelemtwise division by M.
     *
     * @param  M        The divisor.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>
    operator/          (const Matrix<S,P>& M) const {
		Matrix<T,P> res = *this;
		return res /= M;
    }


    /**
     * @brief           Elementwise division by scalar. i.e. this * m.
     *
     * @param  s        The divisor.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>
    operator/           (const S s) const {
		Matrix<T,P> res = *this;
		return res /= s;
	}


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator/=         (const Matrix<T,P>& M) {
        op_assert (_dim == M.Dim(), *this, M);
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::divides<T>());
        return *this;
    }


    /**
     * @brief           ELementwise division and assignment operator. i.e. this = this ./ M.
     *
     * @param  M        Divisor matrix.
     * @return          Result
     */
    template <class S> inline Matrix<T,P>&
    operator/=         (const Matrix<S,P> &M) {
        op_assert (_dim == M.Dim(), *this, M);
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::divides<T>());
        return *this;
    }



    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Divisor scalar.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>&
    operator/=         (const S s) {
        std::transform (_M.begin(), _M.end(), _M.begin(), std::bind2nd(std::divides<T>(),(T)s));
        return *this;
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
    operator*  (const double s, const Matrix<T,P>& m) {
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
    operator*  (const float s, const Matrix<T,P> &m) {
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
    operator*  (const short s, const Matrix<T,P> &m) {
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
    operator*  (const int s, const Matrix<T,P> &m) {
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
    operator*  (const long s, const Matrix<T,P> &m) {
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
    operator*  (const cxfl s, const Matrix<T,P> &m) {
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
    operator*  (const cxdb s, const Matrix<T,P> &m) {
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
    operator+  (const double s, const Matrix<T,P> &m) {
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
    operator+  (const float s, const Matrix<T,P> &m) {
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
    operator+  (const short s, const Matrix<T,P> &m) {
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
    operator+  (const int s, const Matrix<T,P> &m) {
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
    operator+  (const long s, const Matrix<T,P> &m) {
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
    operator+  (const cxfl s, const Matrix<T,P> &m) {
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
    operator+  (const cxdb s, const Matrix<T,P> &m) {
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
    operator-  (const double s, const Matrix<T,P> &m) {
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
    operator-  (const float s, const Matrix<T,P> &m) {
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
    operator-  (const short s, const Matrix<T,P> &m) {
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
    operator-  (const int s, const Matrix<T,P> &m) {
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
    operator-  (const long s, const Matrix<T,P> &m) {
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
    operator-  (const cxfl s, const Matrix<T,P> &m) {
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
    operator-  (const cxdb s, const Matrix<T,P> &m) {
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
    operator/  (const T s, const Matrix<T,P> &m) {
        Matrix<T,P> res = m;
        for (size_t i = 0; i < m.Size(); ++i)
            res[i] = s / res[i];
        return res;
    }

    /**
     * @brief           Elementwise equality with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m == s
     */
    inline friend Matrix<cbool>
    operator== (const T s, const Matrix<T,P>& m) {
        return   m == s;
    }


    /**
     * @brief           Elementwise >= with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m <= t
     */
    inline friend Matrix<cbool>
    operator>= (const T s, const Matrix<T,P>& m) {
        return   m <= s;
    }


    /**
     * @brief           Elementwise <= with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T<=M
     */
    inline friend Matrix<cbool>
    operator<= (const T s, const Matrix<T,P>& m) {
        return   m >= s;
    }


    /**
     * @brief           Elementwise unequality with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T!=M
     */
    inline friend Matrix<cbool>
    operator!= (const T s, const Matrix<T,P>& m) {
        return   m != s;
    }


    /**
     * @brief           Elementwise equality with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T+M
     */
    inline friend Matrix<cbool>
    operator>  (const T s, const Matrix<T,P>& m) {
        return   m <  s;
    }


    /**
     * @brief           Elementwise < with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T+M
     */
    inline friend Matrix<cbool>
    operator<  (const T s, const Matrix<T,P>& m) {
        return   m >  s;
    }


    /**
     * @brief           Elementwise equality with scalar (lhs)
     *
     * @param  mb       Scalar lhs
     * @param  m        Matrix rhs
     * @return          T+M
     */
    //inline friend Matrix<T,P>
    //operator&  (const Matrix<cbool>& mb, const Matrix<T,P>& m) {
    //    return   m & mb;
// }

    //@}


    /**
     * @brief           Elementwise equality, result[i] = (this[i] == m[i]). i.e. this == m
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator==          (const Matrix<T,P>& M) const {
        op_assert (_dim == M.Dim(), *this, M);
        Matrix<cbool> res(_dim);
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
        op_assert (_dim == M.Dim(), *this, M);
		Matrix<cbool> res (_dim);
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
		for (size_t i = 0; i < Size(); ++i)
			res[i] =  (_M[i] == s) ? 1 : 0;
        return res;
    }



protected:
	
    /**
     * @brief           Number of elements
     *
     * @return          Size
     */
    inline size_t
    DimProd () const {
      return std::accumulate(_dim.begin(), _dim.end(), size_t(1), c_multiply<size_t>);
    }

    /**
     * @brief          Allocate RAM
     */
    inline void
    Allocate () {
        size_t ds = _dim.size(), i;
		_dsz.resize(ds,1);
	    for (i = 1; i < ds; ++i)
	        _dsz[i] = _dsz[i-1]*_dim[i-1];
        size_t  n = DimProd(); 
        if (n != _M.size())
            _M.resize(n);
    }

    // Structure
    Vector<size_t> _dim; /// Dimensions
    Vector<size_t> _dsz; /// Dimension size.
    Vector<float>  _res; /// Resolutions
    
    //Data
    Vector<T>        _M; /// Data container

    // Name
    std::string         _name; /// Name
    
#ifdef HAVE_SCALAPACK
    // BLACS 
	int              _bs;
	int              _desc[9]; /**< @brief matrix grid vector */
	int              _gdim[2]; /**< @brief Global dimensions */
#endif
    
};

#endif // __MATRIX_H__


