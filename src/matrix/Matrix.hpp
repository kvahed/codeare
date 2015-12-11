/*
 *  codeare Copyright (C) 2010-2014
 *                        Kaveh Vahedipour
 *                        Forschungszentrum Juelich, Germany
 *                        NYU School of Medicine, New York, USA
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
#include "SIMDTraits.hpp"

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
#include <typeinfo>

#ifdef HAS_CXX11_TUPLE
	#include <tuple>
#else
	#include <boost/tuple/tuple.hpp>
#endif

#include <boost/lexical_cast.hpp>

#include <Assert.hpp>

static inline std::vector<std::string> Parse (const std::string& str, const std::string& dlm) {
	assert (dlm.size() > 0);
	std::vector<std::string> sv;
	size_t  start = 0, end = 0;
	while (end != std::string::npos) {
		end = str.find (dlm, start);
		sv.push_back(str.substr(start, (end == std::string::npos) ?
				std::string::npos : end - start));
		start = ((end > (std::string::npos - dlm.size())) ?
				std::string::npos : end + dlm.size());
	}
	return sv;
}

static const int end = -1;

// Pretty print function names
#if (0 < _MSC_VER)
  #define PRETTY_FUNCTION __FUNCSIG__
#else
  #define PRETTY_FUNCTION __PRETTY_FUNCTION__
#endif


#include "View.hpp"

enum MatrixException {
    DIMS_VECTOR_EMPTY, //0
    DIMS_VECTOR_CONTAINS_ZEROS, //1
	MUST_HAVE_MATCHING_DIMENSIONS_AND_RESOLUTIONS_VECTORS, //2
	ZERO_SIDE_LENGTH, //3
	ZERO_NUMBER_COLUMNS, //4
	ZERO_NUMBER_ROWS, //5
	ZERO_NUMBER_SLICES, //6
	INDEX_EXCEEDS_NUMBER_ELEMENTS, //7
    DIMENSIONS_MUST_MATCH, //8
	INDEX_EXCEEDS_DIMENSION, //9
	DIMENSION_ECXEEDS_DIMENSIONALITY, //10
	CONTAINER_SIZE_MUST_MATCH, //11
	TWO_DIMENSIONAL_OPERATION, //12
	NEGATIVE_INDEX, //13
    REDIM_MUST_PRESERVE_SIZE
};

static const char* MatrixExceptionMessages[] = {
    "Empty dimensions vector", //0
    "Dimensions vector contains 0s" //1
	"Specified dimensions and resolutions vectors have different lengths", //2
	"Matrix with zero side length", //3
	"Matrix with zero height", //4
	"Matrix with zero width", //5
	"Matrix with zero slices", //6
	"Index exceeds number of elements", //7
	"Dimensions must match", //8
	"Index exceeds dimension", //9
	"Dimension exceeds dimensionality", //10
	"Container size must match", //11
	"2D operation only", //12
	"Negative index", //13
    "Redimensioning must preserve size" //14
};

inline static void report_and_throw (const char* fname, const size_t& lnumber,
                                     const char*  func, const MatrixException& x, 
                                     const long&   n1 = -1, const long& n2 = -1) {
    std::cerr << fname << ":" << lnumber << "\n \t" << func << "\n \t"
              << "*** ERROR: " << MatrixExceptionMessages[x-1];
    if (n1 > -1)
        std::cerr << " n1: " << n1;
    if (n2 > -1)
        std::cerr << " n2: " << n2;
	std::cerr << std::endl;
    throw x;
}

#ifndef DNDEBUG
#  ifndef MATRIX_ASSERT
#    define MATRIX_ASSERT(c,x) if (!(c))				\
		report_and_throw (__FILE__, __LINE__, PRETTY_FUNCTION, x)
#  endif
#  ifndef MATRIX_ASSERT2
#    define MATRIX_ASSERT2(c,x,n1,n2) if (!(c))				\
       report_and_throw (__FILE__, __LINE__, PRETTY_FUNCTION, x, n1, n2)
#  endif
#endif

#ifdef EXPLICIT_SIMD
#    include "SIMD.hpp"
#endif

/**
 * @brief   Matrix template.<br/>
 *          Core data structure
 *
 * @author  Kaveh Vahedipour
 * @date    Mar 2010
 */
template <class T, paradigm P> class Matrix : public MatrixType<T,P> {

public:

    typedef View<T,true>  RHSView;
    typedef View<T,false> LHSView;
    
    /**
     * @name Constructors and destructors
     *       Constructors and destructors
     */
    //@{
    
    
    /**
     * @brief           Contruct 1-dim with single element.
     */
	inline Matrix () NOEXCEPT {
        _dim.resize(1,1);
        _res.resize(1,1.0);
        Allocate();
    }
	
	
    /**
     * @brief           Construct matrix with aligned dimension vector
     *
     * @param  dim      All dimensions
     */
	inline Matrix (const Vector<size_t>& dim) {
	    _dim = dim;
		MATRIX_ASSERT(!_dim.empty(), DIMS_VECTOR_EMPTY);
        MATRIX_ASSERT(std::find(dim.begin(),dim.end(),size_t(0))==dim.end(),
        		DIMS_VECTOR_CONTAINS_ZEROS);
        _res.resize(_dim.size(),1.0);
        Allocate();
	}
	
	
    /**
     * @brief           Construct matrix with aligned dimension vector
     *
     * @param  dim      All dimensions
     */
	inline Matrix (const size_t* dims, const size_t ndims) {
	    _dim.resize(ndims);
	    std::copy(dims, dims+ndims, _dim.begin());
		MATRIX_ASSERT(!_dim.empty(), DIMS_VECTOR_EMPTY);
        MATRIX_ASSERT(std::find(_dim.begin(),_dim.end(),size_t(0))==_dim.end(),
        		DIMS_VECTOR_CONTAINS_ZEROS);
        _res.resize(_dim.size(),1.0);
        Allocate();
	}


    /**
     * @brief           Construct matrix with dimension and resolution arrays
     *
     * @param  dim      All 16 Dimensions
     * @param  res      All 16 Resolutions
     */
	inline explicit Matrix (const Vector<size_t>& dim, const Vector<float>& res) {
		_dim = dim;
        MATRIX_ASSERT(!_dim.empty(),DIMS_VECTOR_EMPTY);
        MATRIX_ASSERT(std::find(dim.begin(),dim.end(),size_t(0))==dim.end(),
            DIMS_VECTOR_CONTAINS_ZEROS);
        MATRIX_ASSERT(dim.size()==res.size(),
        		MUST_HAVE_MATCHING_DIMENSIONS_AND_RESOLUTIONS_VECTORS);
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
    inline explicit Matrix (const size_t& n) {
    	MATRIX_ASSERT(n!=0,ZERO_SIDE_LENGTH);
		_dim.resize(2,n);
    	std::cout << _dim << std::endl;
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
    inline Matrix (const size_t& m, const size_t& n)  {
    	MATRIX_ASSERT(n!=0,ZERO_NUMBER_COLUMNS);
    	MATRIX_ASSERT(m!=0,ZERO_NUMBER_ROWS);
    	_dim.resize(2); _dim[0] = m; _dim[1] = n;
    	std::cout << _dim << std::endl;
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
    inline Matrix (const size_t& m, const size_t& n, const size_t& k)  {
   		MATRIX_ASSERT(n!=0,ZERO_NUMBER_COLUMNS);
   		MATRIX_ASSERT(m!=0,ZERO_NUMBER_ROWS);
   		MATRIX_ASSERT(k!=0,ZERO_NUMBER_SLICES);
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
    inline Matrix (const size_t& col,     const size_t& lin,     const size_t& cha,     const size_t& set,
                   const size_t& eco = 1, const size_t& phs = 1, const size_t& rep = 1, const size_t& seg = 1,
                   const size_t& par = 1, const size_t& slc = 1, const size_t& ida = 1, const size_t& idb = 1,
                   const size_t& idc = 1, const size_t& idd = 1, const size_t& ide = 1, const size_t& ave = 1)  {

		MATRIX_ASSERT(col && lin && cha && set && eco && phs && rep && seg &&
				par && slc && ida && idb && idc && idd && ide && ave, ZERO_SIDE_LENGTH);

	    size_t nd = 16, n = nd, i;
	    _dim.resize(n);
	    _dim[ 0] = col; _dim[ 1] = lin; _dim[ 2] = cha; _dim[ 3] = set;
	    _dim[ 4] = eco; _dim[ 5] = phs; _dim[ 6] = rep; _dim[ 7] = seg;
	    _dim[ 8] = par; _dim[ 9] = slc; _dim[10] = ida; _dim[11] = idb;
	    _dim[12] = idc; _dim[13] = idd; _dim[14] = ide; _dim[15] = ave;

        // Remove trailing singleton dimensions except of first :)
        for (i = 1; i < nd; i++)
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
    inline Matrix (const Matrix<T,P> &M) NOEXCEPT {
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
    inline Matrix (Matrix<T,P>&& M) NOEXCEPT {
    	if (this != &M)
    		*this = M;
    }

    inline virtual ~Matrix() {}
#endif

#ifdef HAVE_CXX11_CONDITIONAL
    inline Matrix (RHSView& v) {
		_dim = v._dim;
        MATRIX_ASSERT(!_dim.empty(),DIMS_VECTOR_EMPTY);
        MATRIX_ASSERT(std::find(_dim.begin(),_dim.end(),size_t(0))==_dim.end(),
        		DIMS_VECTOR_CONTAINS_ZEROS);
        Allocate();
        for (size_t i = 0; i < Size(); ++i)
            _M[0] = *(v._pointers[i]);
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
    inline const T& operator[] (const size_t& p) const {
        MATRIX_ASSERT2(p<Size(),INDEX_EXCEEDS_NUMBER_ELEMENTS,p,Size());
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
    inline T& operator[] (const size_t& p) {
        MATRIX_ASSERT2(p<Size(),INDEX_EXCEEDS_NUMBER_ELEMENTS,p,Size());
        return _M[p];
    }


    /**
     * @brief           Get pointer to memory starting at p-th (default:0) element.
     *  
     * @param  p        Position
     *
     * @return          Data 
     */
    inline const T* Ptr (const size_t& p = 0) const {
        MATRIX_ASSERT(p<Size(),INDEX_EXCEEDS_NUMBER_ELEMENTS);
        return _M.ptr(p);
    }

    /**
     * @brief           Get pointer to memory starting at p-th (default:0) element.
     *
     * @param  p        Position
     *
     * @return          Data
     */
    inline T* Ptr (const size_t& p = 0) {
        MATRIX_ASSERT(p<Size(),INDEX_EXCEEDS_NUMBER_ELEMENTS);
        return _M.ptr(p);
    }


    
    /**
     * @brief           Data container (lhs)
     *  
     * @return          Data container
     */
    inline Vector<T>& Container () {
        return _M;
    }

    
    /**
     * @brief           Data container (rhs)
     *  
     * @return          Data container
     */
    inline Vector<T> Container () const {
        return _M;
    }


    /**
     * @brief           Container iterator to first element (lhs)
     *
     * @return          Container iterator
     */
    inline typename Vector<T>::iterator Begin () {
    	return _M.begin ();
    }


    /**
     * @brief           Container const iterator to first element (rhs)
     *
     * @return          Container const iterator
     */
    inline typename Vector<T>::const_iterator Begin () const {
    	return _M.begin ();
    }

    
    /**
     * @brief           Container iterator to last element (lhs)
     *
     * @return          Container iterator
     */
    inline typename Vector<T>::iterator End () {
    	return _M.end ();
    }


    /**
     * @brief           Container const iterator to last element (rhs)
     *
     * @return          Container const iterator
     */
    inline typename Vector<T>::const_iterator End () const {
    	return _M.end ();
    }


    /**
     * @brief           Element at position p (rhs)
     *  
     * @param  p        Position
     * @return          Value at _M[p]
     */
    inline const T& At (const size_t& p) const {
        MATRIX_ASSERT(p<Size(),INDEX_EXCEEDS_NUMBER_ELEMENTS);
        return _M[p];
    }


    /**
     * @brief            Element at position (lhs)
     *  
     * @param  pos       Position
     * @return           Reference to _M[p]
     */
    inline T& At (const size_t& p) {
        MATRIX_ASSERT(p<Size(),INDEX_EXCEEDS_NUMBER_ELEMENTS);
        return _M[p];
    }


    
    /**
     * @brief           Get element in (first) slice (rhs)
     *  
     * @param  x        Column
     * @param  y        Line
     *
     * @return          Value
     */
    inline const T& At (const size_t& x, const size_t& y) const {
        MATRIX_ASSERT2(x<_dim[0],INDEX_EXCEEDS_DIMENSION,x,_dim[0]);
        MATRIX_ASSERT2(y==0 || y<_dim[1],INDEX_EXCEEDS_DIMENSION,y,_dim[1]);
        return _M[x + _dim[0]*y];
    }
    inline       T& At (const size_t& x, const size_t& y)       {
        MATRIX_ASSERT2(x<_dim[0],INDEX_EXCEEDS_DIMENSION,x,_dim[0]);
        MATRIX_ASSERT2(y==0 || y<_dim[1],INDEX_EXCEEDS_DIMENSION,y,_dim[1]);
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
    inline const T& At (const size_t& x, const size_t& y, const size_t& z) const {
        MATRIX_ASSERT(x<_dim[0],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(y<_dim[1],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(z<_dim[2],INDEX_EXCEEDS_DIMENSION);
        return _M[x + _dsz[1]*y + _dsz[2]*z];
    }
    inline T&       At (const size_t& x, const size_t& y, const size_t& z)  {
        MATRIX_ASSERT(x<_dim[0],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(y<_dim[1],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(z<_dim[2],INDEX_EXCEEDS_DIMENSION);
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
    inline const T& At (const size_t& x, const size_t& y, const size_t& z,
    		const size_t& w) const  {
        MATRIX_ASSERT(x<_dim[0],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(y<_dim[1],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(z<_dim[2],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(w<_dim[3],INDEX_EXCEEDS_DIMENSION);
        return _M[x + _dsz[1]*y + _dsz[2]*z + _dsz[3]*w];
    }
    inline T& At (const size_t& x, const size_t& y, const size_t& z,
    		const size_t& w)  {
        MATRIX_ASSERT(x<_dim[0],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(y<_dim[1],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(z<_dim[2],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(w<_dim[3],INDEX_EXCEEDS_DIMENSION);
        return _M[x + _dsz[1]*y + _dsz[2]*z + _dsz[3]*w];
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
    inline const T& At (const size_t& n0, const size_t& n1, const size_t& n2,
                        const size_t& n3, const size_t& n4) const  {
        MATRIX_ASSERT(n0<_dim[0],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n1<_dim[1],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n2<_dim[2],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n3<_dim[3],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n4<_dim[4],INDEX_EXCEEDS_DIMENSION);
        return _M[n0 + _dsz[1]*n1 + _dsz[2]*n2 + _dsz[3]*n3 + _dsz[4]*n4];
    }
    inline       T& At (const size_t& n0, const size_t& n1, const size_t& n2,
                        const size_t& n3, const size_t& n4) {
        MATRIX_ASSERT(n0<_dim[0],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n1<_dim[1],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n2<_dim[2],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n3<_dim[3],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n4<_dim[4],INDEX_EXCEEDS_DIMENSION);
        return _M[n0 + _dsz[1]*n1 + _dsz[2]*n2 + _dsz[3]*n3 + _dsz[4]*n4];
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
    inline const T& At (const size_t& n0, const size_t& n1, const size_t& n2,
                        const size_t& n3, const size_t& n4, const size_t& n5) const  {
        MATRIX_ASSERT(n0<_dim[0],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n1<_dim[1],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n2<_dim[2],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n3<_dim[3],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n4<_dim[4],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n5<_dim[5],INDEX_EXCEEDS_DIMENSION);
        return _M[n0 + _dsz[1]*n1 + _dsz[2]*n2 + _dsz[3]*n3 + _dsz[4]*n4 + _dsz[5]*n5];
    }
    inline       T& At (const size_t& n0, const size_t& n1, const size_t& n2,
                        const size_t& n3, const size_t& n4, const size_t& n5) {
        MATRIX_ASSERT(n0<_dim[0],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n1<_dim[1],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n2<_dim[2],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n3<_dim[3],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n4<_dim[4],INDEX_EXCEEDS_DIMENSION);
        MATRIX_ASSERT(n5<_dim[5],INDEX_EXCEEDS_DIMENSION);
        return _M[n0 + _dsz[1]*n1 + _dsz[2]*n2 + _dsz[3]*n3 + _dsz[4]*n4+ _dsz[5]*n5];
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
    inline const T& At
		(const size_t& n00,     const size_t& n01,	   const size_t& n02,
		 const size_t& n03,     const size_t& n04,     const size_t& n05,
		 const size_t& n06,     const size_t& n07 = 0, const size_t& n08 = 0,
		 const size_t& n09 = 0, const size_t& n10 = 0, const size_t& n11 = 0,
		 const size_t& n12 = 0, const size_t& n13 = 0, const size_t& n14 = 0,
		 const size_t& n15 = 0) const  {
    	MATRIX_ASSERT (n00<_dim[ 0] && n01<_dim[ 1] && n02<_dim[ 2] && n03<_dim[ 3]
				    && n04<_dim[ 4] && n05<_dim[ 5] && n06<_dim[ 6] && n07<_dim[ 7]
					&& n08<_dim[ 8] && n09<_dim[ 9] && n10<_dim[10] && n11<_dim[11]
                    && n12<_dim[12] && n13<_dim[13] && n14<_dim[14] && n15<_dim[15],
					   INDEX_EXCEEDS_DIMENSION);
        return _M [n00          + n01*_dsz[ 1] + n02*_dsz[ 2] + n03*_dsz[ 3] +
				   n04*_dsz[ 4] + n05*_dsz[ 5] + n06*_dsz[ 6] + n07*_dsz[ 7] +
				   n08*_dsz[ 8] + n09*_dsz[ 9] + n10*_dsz[10] + n11*_dsz[11] +
				   n12*_dsz[12] + n13*_dsz[13] + n14*_dsz[14] + n15*_dsz[15]];
    }
    inline T& At
		(const size_t& n00,     const size_t& n01,	   const size_t& n02,
   		 const size_t& n03,     const size_t& n04,     const size_t& n05,
   		 const size_t& n06,     const size_t& n07 = 0, const size_t& n08 = 0,
   		 const size_t& n09 = 0, const size_t& n10 = 0, const size_t& n11 = 0,
   		 const size_t& n12 = 0, const size_t& n13 = 0, const size_t& n14 = 0,
   		 const size_t& n15 = 0) {
    	MATRIX_ASSERT (n00<_dim[ 0] && n01<_dim[ 1] && n02<_dim[ 2] && n03<_dim[ 3]
				    && n04<_dim[ 4] && n05<_dim[ 5] && n06<_dim[ 6] && n07<_dim[ 7]
					&& n08<_dim[ 8] && n09<_dim[ 9] && n10<_dim[10] && n11<_dim[11]
                    && n12<_dim[12] && n13<_dim[13] && n14<_dim[14] && n15<_dim[15],
					   INDEX_EXCEEDS_DIMENSION);
           return _M [n00          + n01*_dsz[ 1] + n02*_dsz[ 2] + n03*_dsz[ 3] +
   				   n04*_dsz[ 4] + n05*_dsz[ 5] + n06*_dsz[ 6] + n07*_dsz[ 7] +
   				   n08*_dsz[ 8] + n09*_dsz[ 9] + n10*_dsz[10] + n11*_dsz[11] +
   				   n12*_dsz[12] + n13*_dsz[13] + n14*_dsz[14] + n15*_dsz[15]];
    }
    

    /**
     * @brief          Cast operator
     *
     * @return         Casted copy
     */
    template<class S> inline operator Matrix<S,P> () const {
		Matrix<S,P> m(_dim);
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
    inline const T& operator() (const size_t& p) const {
        return this->At(p);
    }

    
    /**
     * @brief           Get value of pth element of repository.
     *
     * @param  p        Requested position.
     * @return          Requested scalar value.
     */
    inline T& operator() (const size_t& p) {
        return this->At(p);
    }

    
    /**
     * @brief           Get value in slice
     *
     * @param  x        Column
     * @param  y        Line
     * @return          Value
     */
    inline const T& operator() (const size_t& x, const size_t& y) const {
        return this->At(x,y);
    }
    

    /**
     * @brief           Reference to value in slice
     *
     * @param  x        Column
     * @param  y        Line
     * @return          Reference
     */
    inline T& operator() (const size_t& x, const size_t& y)  {
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
    inline const T& operator() (const size_t& x, const size_t& y, const size_t& z) const {
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
    inline T& operator() (const size_t& x, const size_t& y, const size_t& z) {
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
    inline const T& operator() (const size_t& x, const size_t& y, const size_t& z,
                                const size_t& w) const {
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
    inline T& operator() (const size_t& x, const size_t& y, const size_t& z,
    		const size_t& w) {
        return this->At(x,y,z,w);
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
    inline const T& operator() (const size_t& n0, const size_t& n1, const size_t& n2,
                                const size_t& n3, const size_t& n4) const {
        return this->At(n0,n1,n2,n3,n4);
    }
    inline T& operator() (const size_t& n0, const size_t& n1, const size_t& n2,
                          const size_t& n3, const size_t& n4) {
        return this->At(n0,n1,n2,n3,n4);
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
    inline const T& operator() (const size_t& n0, const size_t& n1, const size_t& n2,
                                const size_t& n3, const size_t& n4, const size_t& n5) const {
        return this->At(n0,n1,n2,n3,n4,n5);
    }
    inline T& operator() (const size_t& n0, const size_t& n1, const size_t& n2,
                          const size_t& n3, const size_t& n4, const size_t& n5) {
        return this->At(n0,n1,n2,n3,n4,n5);
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
    inline T& operator()
    		(const size_t& n00,     const size_t& n01,	   const size_t& n02,
      		 const size_t& n03,     const size_t& n04,     const size_t& n05,
      		 const size_t& n06,     const size_t& n07 = 0, const size_t& n08 = 0,
      		 const size_t& n09 = 0, const size_t& n10 = 0, const size_t& n11 = 0,
      		 const size_t& n12 = 0, const size_t& n13 = 0, const size_t& n14 = 0,
      		 const size_t& n15 = 0) {
        return _M [n00          + n01*_dsz[ 1] + n02*_dsz[ 2] + n03*_dsz[ 3] +
				   n04*_dsz[ 4] + n05*_dsz[ 5] + n06*_dsz[ 6] + n07*_dsz[ 7] +
				   n08*_dsz[ 8] + n09*_dsz[ 9] + n10*_dsz[10] + n11*_dsz[11] +
				   n12*_dsz[12] + n13*_dsz[13] + n14*_dsz[14] + n15*_dsz[15]];
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
    inline const T& operator()
    		(const size_t& n00,     const size_t& n01,	   const size_t& n02,
     		 const size_t& n03,     const size_t& n04,     const size_t& n05,
     		 const size_t& n06,     const size_t& n07 = 0, const size_t& n08 = 0,
     		 const size_t& n09 = 0, const size_t& n10 = 0, const size_t& n11 = 0,
     		 const size_t& n12 = 0, const size_t& n13 = 0, const size_t& n14 = 0,
     		 const size_t& n15 = 0) const {
       return _M [n00          + n01*_dsz[ 1] + n02*_dsz[ 2] + n03*_dsz[ 3] +
				  n04*_dsz[ 4] + n05*_dsz[ 5] + n06*_dsz[ 6] + n07*_dsz[ 7] +
				  n08*_dsz[ 8] + n09*_dsz[ 9] + n10*_dsz[10] + n11*_dsz[11] +
				  n12*_dsz[12] + n13*_dsz[13] + n14*_dsz[14] + n15*_dsz[15]];
	}


    //@}
    

    
    
#ifndef NO_LAPACK
    /**
     * @brief           Matrix product. i.e. this * M.
     *
     * @param  M        The factor.
     */
    Matrix<T,P> operator->* (const Matrix<T,P>& M) const;

    /**
     * @brief           Matrix Product.
     *
     * @param   M       The factor
     * @param   transa  Transpose ('T') / Conjugate transpose ('C') the left matrix. Default: No transposition'N'
     * @param   transb  Transpose ('T') / Conjugate transpose ('C') the right matrix. Default: No transposition 'N'
     * @return          Product of this and M.
     */
    Matrix<T,P> mul (const Matrix<T,P> &M, const char transa = 'N', const char transb = 'N') const;



    /**
     * @brief           Complex conjugate left and multiply with right.
     *
     * @param   M       Factor
     * @return          Product of conj(this) and M.
     */
    Matrix<T,P> mult  (const Matrix<T,P> &M) const;


    /**
     * @brief           Scalar product (complex: conjugate first vector) using
     *                  <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTC and XDOT
     *
     * @param  M        Factor
     * @return          Scalar product
     */
     T dotc (const Matrix<T,P>& M) const;


    /**
     * @brief           Scalar product using <a href="http://www.netlib.org/blas/">BLAS</a> routines XDOTU and XDOT
     *
     * @param  M        Factor
     * @return          Scalar product
     */
    T dot (const Matrix<T,P>& M) const;
    
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
    inline size_t Height () const NOEXCEPT {
    	return _dim[0];
    }


    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    inline size_t Width () const NOEXCEPT {
        return (_dim.size()>1) ? _dim[1] : 1;
    }

#ifdef HAVE_SCALAPACK

    /**
     * @brief           Get number of rows, i.e. tmp = size(this); tmp(1).
     *
     * @return          Number of rows.
     */
    inline size_t GHeight () const {
        return _gdim[0];
    }


    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    inline size_t GWidth () const {
        return _gdim[1];
    }


    /**
     * @brief           Get number of columns, i.e. tmp = size(this); tmp(2).
     *
     * @return          Number of columns.
     */
    inline const int* Desc () const NOEXCEPT {
        return _desc;
    }
#endif


    /**
     * @brief           Get resolution of a given dimension.
     *
     * @param   i       Dimension
     * @return          Resolution .
     */
    inline float Res (const size_t& i) const {
        MATRIX_ASSERT (i < _dim.size(), DIMENSION_ECXEEDS_DIMENSIONALITY);
        return _res[i];
    }


    /**
     * @brief           Resolution of a given dimension.
     *
     * @param   i       Dimension
     * @return          Resolution
     */
    inline float& Res (const size_t& i) {
        MATRIX_ASSERT (i < _dim.size(), DIMENSION_ECXEEDS_DIMENSIONALITY);
        return _res[i];
    }



    /**
     * @brief           Resolution array
     *
     * @return          All resolutions
     */
    inline const Vector<float>& Res () const NOEXCEPT {
        return _res;
    }



    /**
     * @brief           Get size a given dimension.
     *
     * @param   i       Dimension
     * @return          Dimension
     */
    inline virtual size_t Dim (const size_t& i) const {
        return (i < _dim.size()) ? _dim[i]: 1;
    }


    /**
     * @brief           Get dimension array
     *
     * @return          All dimensions
     */
    inline virtual const Vector<size_t>& Dim() const NOEXCEPT {
        return _dim;
    }


    inline void Squeeze () {
        for (auto i = _dim.begin(), j = _dsz.begin(); i != _dim.end(); ++i, ++j)
            if (*i == 1) {
                _dsz.erase(j--);
                _dim.erase(i--);
            }
    }



    /**
     * @brief           Number of dimensions
     * @return          Number of dimensions
     */
    inline size_t NDim() const NOEXCEPT {
    	return _dim.size();
    }


    /**
     * @brief           Get dimension sizes
     *
     * @return          All dimensions
     */
    inline const Vector<size_t>& Dsz () const NOEXCEPT {
        return _dsz;
    }


    /**
     * @brief           Purge data and free RAM.
     */
    inline void Clear() NOEXCEPT {
    	_dim.clear();
        _dsz.clear();
        _res.clear();
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
    inline const char* GetClassName () const NOEXCEPT {
        return _name.c_str(); 
    }

    /**
     * @brief           Who are we?
     *
     * @return          Class name
     */ 
    inline void SetClassName (const char* name) NOEXCEPT {
        _name = name; 
    }


    /**
     * @brief           Number of elements
     *
     * @return          Size
     */
    inline virtual size_t Size () const NOEXCEPT {
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
    inline Matrix<T,P>& operator= (Matrix<T,P>&& rhs) NOEXCEPT {
        _M    = std::move(rhs._M);
        _name = std::move(rhs._name);
        _res  = std::move(rhs._res);
        _dsz  = std::move(rhs._dsz);
        _dim  = std::move(rhs._dim);
        return *this;
    }
	/**
	 * @brief  Copy assignment operator
	 * @param  rhs  Right hand side
	 * @return      Left hand side
	 */
    inline Matrix<T,P>& operator= (const Matrix<T,P>& rhs) NOEXCEPT {
        _M    = rhs._M;
        _name = rhs._name;
        _res  = rhs._res;
        _dsz  = rhs._dsz;
        _dim  = rhs._dim;
        return *this;
    }
#endif
    /**
     * @brief           Assignment data
     *
     * @param  v        Data vector (size must match numel(M)).
     */
    inline Matrix<T,P>& operator= (const Vector<T>& v) {
    	if (_M.size() == 1) { // we are being assigned out of nothing
    		_dim.resize(1,v.size());
    	    _res.resize(1,1.0);
    	    Allocate();
    	} else {              // we have been allocated already
    		MATRIX_ASSERT(_M.size() == v.size(), CONTAINER_SIZE_MUST_MATCH);
    	}
    	if (&_M != &v)
            _M = v;
        return *this;
    }

#ifdef HAVE_CXX11_CONDITIONAL
    template<class S>
    inline Matrix<T,P>& operator= (const View<S,true>& v) {
        _dim = v._dim;
        Allocate();
        for (size_t i = 0; i < Size(); ++i)
            _M[i] = v[i];
        return *this;
    }
#endif


    /**
     * @brief           Assignment operator. Sets all elements s.
     *
     * @param  s        The assigned scalar.
     */
    inline const Matrix<T,P>& operator= (const T& s) {
        std::fill (_M.begin(), _M.end(), s);
        return *this;
    }


    /**
     * @brief           Unary minus (additive inverse)
     *
     * @return          Negation
     */
    inline Matrix<T,P> operator- () const {
        Matrix<T,P> res (_dim);
        std::transform (_M.begin(), _M.end(), res.Begin(), std::negate<T>());
        return res;
    }


    /**
     * @brief           Unary plus
     *
     * @return          Identity
     */
    inline Matrix<T,P> operator+ () const {
        return *this;
    }


    /**
     * @brief           Complex conjugation. i.e. this.'
     *
     * @return          Matrix::tr()
     */
    inline Matrix<T,P> operator! () const {
    	Matrix<T,P> res = *this;
    	Vec(_M, res._M, codeare::conjugate<T>());
    	return res;
    }


     /**
     * @brief           Scalar inequality. result[i] = (this[i] != m). i.e. this ~= m
     *
     * @param  s        Comparing scalar.
     * @return          Matrix of false where elements are equal s and true else.
     */
    inline Matrix<cbool> operator!= (const T& s) const {
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
    inline Matrix<cbool> operator> (const T& s) const {
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
    inline Matrix<cbool> operator>= (const T& s) const {
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
    inline Matrix<cbool> operator<= (const T& s) const {
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
    inline Matrix<cbool> operator< (const T& s) const {
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
    template<class S>
    inline Matrix<cbool> operator!= (const Matrix<S,P>& M) const {
        MATRIX_ASSERT (_dim == M.Dim(), DIMENSIONS_MUST_MATCH);
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
    template<class S>
    inline Matrix<cbool>operator>= (const Matrix<S,P>& M) const {
    	MATRIX_ASSERT (_dim == M.Dim(), DIMENSIONS_MUST_MATCH);
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
    template<class S>
    inline Matrix<cbool> operator<= (const Matrix<S,P>& M) const {
    	MATRIX_ASSERT (_dim == M.Dim(), DIMENSIONS_MUST_MATCH);
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
    template<class S>
    inline Matrix<cbool> operator> (const Matrix<S,P>& M) const {
    	MATRIX_ASSERT (_dim == M.Dim(), DIMENSIONS_MUST_MATCH);
        Matrix<cbool> res(_dim);
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
    template<class S>
    inline Matrix<cbool> operator< (const Matrix<S,P>& M) const {
    	MATRIX_ASSERT (_dim == M.Dim(), DIMENSIONS_MUST_MATCH);
        Matrix<cbool> res(_dim);
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
    template<class S>
    inline Matrix<cbool> operator|| (const Matrix<S,P>& M) const {
    	MATRIX_ASSERT (_dim == M.Dim(), DIMENSIONS_MUST_MATCH);
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
    template<class S>
    inline Matrix<short> operator| (const Matrix<S,P>& rhs) const {
    	MATRIX_ASSERT (_dim == rhs.Dim(), DIMENSIONS_MUST_MATCH);
        Matrix<cbool> ret(_dim);
        for (size_t i = 0; i < Size(); ++i)
        	ret[i] = CompTraits<T>::logical_or(_M[i], ret[i]);
        return ret;
    }


    /**
     * @brief           Matrix comparison, result[i] = (m[i] && this[i] ? 1 : 0). i.e. this & m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    template<class S>
    inline Matrix<cbool> operator&& (const Matrix<S,P>& M) const {
    	MATRIX_ASSERT (_dim == M.Dim(), DIMENSIONS_MUST_MATCH);
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
    template<class S>
    inline Matrix<T,P> operator^ (const S& p) const  {
    	Matrix<T,P> res = *this;
        return res ^= p;
    }


    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Power.
     * @return          Result
     */
    template<class S>
    inline Matrix<T,P>& operator^= (const S& p) {
		for (size_t i = 0; i < Size(); ++i)
			_M[i] = TypeTraits<T>::Pow(_M[i], p);
        return *this;
    }

    /**
     * @brief           Elementwise multiplication (calls one of the below 3)
     * @param  s        Factor
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P> operator+ (const S& s) const {
        Matrix<T,P> res = *this;
        return res += s;
    }

    /**
     * @brief           Elementwise multiplication and assignment operator with
     * 					other matrix; i.e. A = A.*B.
     * @param  M        Factor matrix.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>& operator+= (const Matrix<S,P>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::plus<T>());
        return *this;
    }
    inline Matrix<T,P>& operator+= (const Matrix<T,P>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        Vec(_M, M._M, _M, codeare::plus<T>());
        return *this;
    }

    /**
     * @brief           Elementwise multiplication and assignment operator with.
     *                  view; i.e. A = A.*B(...,...).
     * @param  M        Factor view.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>& operator+= (const View<S,true>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        for (size_t i = 0; i < Size(); ++i)
            _M[i] += M[i];
        return *this;
    }

    /**
     * @brief           ELementwise multiplication with scalar; i.e. A = A.*s;
     * @param  s        Factor scalar.
     * @return          Result
     */
    inline Matrix<T,P>& operator+= (const T& t) {
        Vec(_M, t, _M, codeare::plus<T>());
		return *this;
    }


    /**
     * @brief           Elementwise multiplication (calls one of the below 3)
     * @param  s        Factor
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P> operator- (const S& s) const {
        Matrix<T,P> res = *this;
        return res -= s;
    }

    /**
     * @brief           Elementwise multiplication and assignment operator with
     * 					other matrix; i.e. A = A.*B.
     * @param  M        Factor matrix.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>& operator-= (const Matrix<S,P>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::minus<T>());
        return *this;
    }
    inline Matrix<T,P>& operator-= (const Matrix<T,P>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        Vec(_M, M._M, _M, codeare::minus<T>());
        return *this;
    }

    /**
     * @brief           Elementwise multiplication and assignment operator with.
     *                  view; i.e. A = A.*B(...,...).
     * @param  M        Factor view.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>& operator-= (const View<S,true>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        for (size_t i = 0; i < Size(); ++i)
            _M[i] -= M[i];
        return *this;
    }

    /**
     * @brief           ELementwise multiplication with scalar; i.e. A = A.*s;
     * @param  s        Factor scalar.
     * @return          Result
     */
    inline Matrix<T,P>& operator-= (const T& t) {
        Vec(_M, t, _M, codeare::minus<T>());
		return *this;
    }

    /**
     * @brief           Elementwise multiplication (calls one of the below 3)
     * @param  s        Factor
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P> operator* (const S& s) const {
        Matrix<T,P> res = *this;
        return res *= s;
    }

    /**
     * @brief           Elementwise multiplication and assignment operator with
     * 					other matrix; i.e. A = A.*B.
     * @param  M        Factor matrix.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>& operator*= (const Matrix<S,P>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::multiplies<T>());
        return *this;
    }
    inline Matrix<T,P>& operator*= (const Matrix<T,P>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        Vec(_M, M._M, _M, codeare::multiplies<T>());
        return *this;
    }


    /**
     * @brief           Elementwise multiplication and assignment operator with.
     *                  view; i.e. A = A.*B(...,...).
     * @param  M        Factor view.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>& operator*= (const View<S,true>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        for (size_t i = 0; i < Size(); ++i)
            _M[i] *= M[i];
        return *this;
    }

    /**
     * @brief           ELementwise multiplication with scalar; i.e. A = A.*s;
     * @param  s        Factor scalar.
     * @return          Result
     */
    inline Matrix<T,P>& operator*= (const T& t) {
        Vec(_M, t, _M, codeare::multiplies<T>());
		return *this;
    }

    /**
     * @brief           Elementwise multiplication (calls one of the below 3)
     * @param  s        Factor
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P> operator/ (const S& s) const {
        Matrix<T,P> res = *this;
        return res /= s;
    }

    /**
     * @brief           Elementwise multiplication and assignment operator with
     * 					other matrix; i.e. A = A.*B.
     * @param  M        Factor matrix.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>& operator/= (const Matrix<S,P>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::divides<T>());
        return *this;
    }
    inline Matrix<T,P>& operator/= (const Matrix<T,P>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        Vec(_M, M._M, _M, codeare::divides<T>());
        return *this;
    }

    /**
     * @brief           Elementwise multiplication and assignment operator with.
     *                  view; i.e. A = A.*B(...,...).
     * @param  M        Factor view.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>& operator/= (const View<S,true>& M) {
        MATRIX_ASSERT (_dim==M.Dim(), DIMENSIONS_MUST_MATCH);
        for (size_t i = 0; i < Size(); ++i)
            _M[i] /= M[i];
        return *this;
    }

    /**
     * @brief           ELementwise multiplication with scalar; i.e. A = A.*s;
     * @param  s        Factor scalar.
     * @return          Result
     */
    inline Matrix<T,P>& operator/= (const T& t) NOEXCEPT {
        Vec(_M, t, _M, codeare::divides<T>());
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
    inline friend Matrix<T,P> operator* (const T& s, const Matrix<T,P> &m)
        NOEXCEPT { return m*s; }


    //--
    /**
     * @brief           Elementwise addition with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          -m + s
     */
    inline friend Matrix<T,P> operator- (const T& s, const Matrix<T,P>& m)
        NOEXCEPT { return -m+s; }

    //--
    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P> operator+ (const T& s, const Matrix<T,P>& m)
        NOEXCEPT { return m+s; }

    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P> operator/ (const T& s, const Matrix<T,P> &m)
        NOEXCEPT {
        Matrix<T,P> res = m;
        for (size_t i = 0; i < m.Size(); ++i)
            res[i] = s / res[i];
        return res;
    }


    //@}


    /**
	 * @brief           Elementwise equality, result[i] = (this[i] == m[i]). i.e. this == m
	 *
	 * @param  M        Comparing matrix.
	 * @return          Hit list
	 */
    template<class S> inline Matrix<cbool> operator== (const Matrix<S,P>& M) const {
        MATRIX_ASSERT (_dim == M._dim, DIMENSIONS_MUST_MATCH);
		Matrix<cbool> res (_dim);
		for (size_t i = 0; i < Size(); ++i)
			res[i] = (_M[i] == M[i]);
		return res;
     }


    /**
     * @brief           Scalar equality. result[i] = (this[i] == m).
     *
     * @param  s        Comparing scalar.
     * @return          Matrix of true where elements are equal s and false else.
     */
    inline Matrix<cbool> operator== (const T& s) const {
    	T t(s);
        Matrix<cbool> res (_dim);
		for (size_t i = 0; i < Size(); ++i)
			res[i] = (_M[i] == s);
        return res;
    }


#ifdef HAVE_CXX11_CONDITIONAL

    inline LHSView operator() (R r) {
        Vector<R> vr;
        vr.push_back (r);
        return LHSView(this, vr);
    }
    inline RHSView operator() (const CR r) const {
        Vector<CR> vr;
        vr.push_back (r);
        return RHSView(this, vr);
    }
    inline LHSView operator() (const R r, const size_t& n) {
        Vector<R> vr;
        vr.push_back (r);
        vr.push_back (R(n));
        return LHSView(this, vr);
    }
    inline RHSView operator() (const CR r, const size_t& n) const {
        Vector<CR> vr;
        vr.push_back (r);
        vr.push_back (CR(n));
        return RHSView(this, vr);
    }
    inline LHSView operator() (R r0, R r1) {
        Vector<R> vr;
        vr.push_back (r0);
        vr.push_back (r1);
        return LHSView(this, vr);
    }
    inline RHSView operator() (CR r0, CR r1) const {
        Vector<CR> vr;
        vr.push_back (r0);
        vr.push_back (r1);
        return RHSView(this, vr);
    }
    inline LHSView operator() (R r0, R r1, R r2) {
        Vector<R> vr;
        vr.push_back (r0);
        vr.push_back (r1);
        vr.push_back (r2);        
        return LHSView(this, vr);
    }
    inline RHSView operator() (CR r0, CR r1, CR r2) const {
        Vector<CR> vr;
        vr.push_back (r0);
        vr.push_back (r1);
        vr.push_back (r2);
        return RHSView(this, vr);
    }
    inline LHSView operator() (R r0, R r1, R r2, R r3) {
        Vector<R> vr;
        vr.push_back (r0);
        vr.push_back (r1);
        vr.push_back (r2);        
        vr.push_back (r3);        
        return LHSView(this, vr);
    }
    inline RHSView operator() (CR r0, CR r1, CR r2, CR r3) const {
        Vector<CR> vr;
        vr.push_back (r0);
        vr.push_back (r1);
        vr.push_back (r2);
        vr.push_back (r3);
        return RHSView(this, vr);
    }
    inline LHSView operator() (R r0, R r1, R r2, R r3, R r4) {
        Vector<R> vr;
        vr.push_back (r0);
        vr.push_back (r1);
        vr.push_back (r2);        
        vr.push_back (r3);        
        vr.push_back (r4);
        return LHSView(this, vr);
    }
    inline RHSView operator() (CR r0, CR r1, CR r2, CR r3, CR r4) const {
        Vector<CR> vr;
        vr.push_back (r0);
        vr.push_back (r1);
        vr.push_back (r2);
        vr.push_back (r3);
        vr.push_back (r4);
        return RHSView(this, vr);
    }
    inline LHSView operator() (R r0, R r1, R r2, R r3, R r4, R r5) {
        Vector<R> vr;
        vr.push_back (r0);
        vr.push_back (r1);
        vr.push_back (r2);        
        vr.push_back (r3);        
        vr.push_back (r4);
        vr.push_back (r5);
        return LHSView(this, vr);
    }
    inline RHSView operator() (CR r0, CR r1, CR r2, CR r3, CR r4, CR r5) const {
        Vector<CR> vr;
        vr.push_back (r0);
        vr.push_back (r1);
        vr.push_back (r2);
        vr.push_back (r3);
        vr.push_back (r4);
        vr.push_back (r4);
        return RHSView(this, vr);
    }
#endif

protected:
	
    /**
     * @brief          Allocate RAM
     */
    inline void Allocate () {
        size_t ds = _dim.size(), i;
		_dsz.resize(ds);
        _dsz[0] = 1;
	    for (i = 1; i < ds; ++i)
	        _dsz[i] = _dsz[i-1]*_dim[i-1];
        size_t n = _dsz.back()*_dim.back(); 
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
    std::string      _name; /// Name
    
#ifdef HAVE_SCALAPACK
    // BLACS 
	int              _bs;
	int              _desc[9]; /**< @brief matrix grid vector */
	int              _gdim[2]; /**< @brief Global dimensions */
#endif
    
};

#endif // __MATRIX_H__


