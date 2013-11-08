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

#include "config.h"
#include "common.h"
#include "OMP.hpp"
#include "Complex.hpp"
#include "SIMD.hpp"

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
    
    
    Matrix ();
    Matrix (const std::vector<size_t>& dim);
	virtual ~Matrix();
    Matrix (const size_t dim[INVALID_DIM]);
    Matrix (const std::vector<size_t>& dim, const std::vector<float>& res);
    Matrix (const size_t dim[INVALID_DIM], const float res[INVALID_DIM]);
    Matrix (const size_t& n);
    Matrix (const size_t& m, const size_t& n);
    Matrix (const size_t& m, const size_t& n, const size_t& k);
    Matrix (const size_t col, 
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
    Matrix (const Matrix<T,P> &M);
    T operator[] (const size_t& p) const;
    T& operator[] (const size_t& p);
    const T* Memory (const size_t p = 0) const;
    inline std::valarray<T>&            
    Container           ();
    inline std::valarray<T>            
    Container           () const;
    inline T            
    At                  (const size_t& p) const;
    inline T&           
    At                  (const size_t& pos) ;
    inline T            
    At                  (const size_t& x, const size_t& y) const;
    virtual const char* GetClassName() const;


        /**
     * @brief           Assignment operator. i.e. this = m.
     *
     * @param  M        The assigned matrix.
     */
    inline Matrix<T,P>&
    operator=           (const Matrix<T,P>& M) {

        if (this != &M) {


            for (size_t i = 0; i < INVALID_DIM; ++i) {
                _dim[i] = M.Dim()[i];
                _res[i] = M.Res()[i];
                _dsz[i] = M.Dsz()[i];
            }

            _M = M.Container();

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
    operator-           (const S& s) const {

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

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

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

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

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
    operator-=          (const S& s) {

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
    inline Matrix<T,P>
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
    operator+           (const S& s) const {

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

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

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

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

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
    operator+=          (const S& s) {

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

        for (size_t i = 2; i < INVALID_DIM; ++i)
            assert (_dim[i] == 1);

        Matrix<T,P> res (_dim[1],_dim[0]);

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
    operator&           (const Matrix<unsigned short>& M) const ;
    
    
     /**
     * @brief           Scalar inequality. result[i] = (this[i] != m). i.e. this ~= m
     *
     * @param  s        Comparing scalar.
     * @return          Matrix of false where elements are equal s and true else.
     */
    inline Matrix<unsigned short>
    operator!=          (const T s) const {

        Matrix<unsigned short> res(_dim);
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
    inline Matrix<unsigned short>
    operator>           (const T s) const {

        Matrix<unsigned short> res(_dim);

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
    inline Matrix<unsigned short>
    operator>=          (const T s) const {

		Matrix<unsigned short> res(_dim);

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
    inline Matrix<unsigned short>
    operator<=          (const T s) const {

        Matrix<unsigned short> res(_dim);

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
    inline Matrix<unsigned short>
    operator<           (const T s) const {

        Matrix<unsigned short> res(_dim);

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
    inline Matrix<unsigned short>
    operator==          (const Matrix<T,P>& M) const {

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

        Matrix<unsigned short> res(_dim);
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
	inline Matrix<unsigned short>
	operator==          (const Matrix<S,P>& M) const {

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

		Matrix<unsigned short> res (_dim);

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
    inline Matrix<unsigned short>
    operator==          (const T s) const {

    	T t = (T) s;

        Matrix<unsigned short> res (_dim);

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
    inline Matrix<unsigned short>
    operator!=          (const Matrix<T,P>& M) const {

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

        Matrix<unsigned short> res(_dim,_res);
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
    inline Matrix<unsigned short>
    operator>=          (const Matrix<T,P>& M) const {

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

        Matrix<unsigned short> res(_dim);

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
    inline Matrix<unsigned short>
    operator<=          (const Matrix<T,P>& M) const {

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

        Matrix<unsigned short> res(_dim);

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
    inline Matrix<unsigned short>
    operator>           (const Matrix<T,P>& M) const {

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

        Matrix<unsigned short> res(_dim,_res);

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
    inline Matrix<unsigned short>
    operator<           (const Matrix<T,P>& M) const {

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

        Matrix<unsigned short> res(_dim,_res);

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
    inline Matrix<unsigned short>
    operator||          (const Matrix<T,P>& M) const {

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

        Matrix<unsigned short> res(_dim);

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
    inline Matrix<unsigned short>
    operator&&          (const Matrix<T,P>& M) const {

        Matrix<unsigned short> res(_dim);

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
    operator*          (const S& s) const  {

        Matrix<T,P> res = *this;
    	T t = T(s);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] *= t;

        return res;

    }



    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator*=         (const Matrix<T,P>& M) {

        size_t i;

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

#if defined USE_VALARRAY
        _M *= M.Container();        
#elif defined EXPLICIT_SIMD
        if (fp_type(_M[0]))
        	SSE::binary<T>(_M, M.Container(), SSE::mul<T>(), _M);
        else
        	for (size_t i = 0; i < Size(); ++i)
        		_M[i] *= M[i];
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
            _M[i] *= M[i];
#endif

        return *this;

    }


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>&
    operator*=         (const Matrix<S,P>& M) {

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			_M[i] *= M[i];

		return *this;

    }


    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this *= s.
     *
     * @param  s        Factor scalar.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>&
    operator*=         (const S& s) {

    	T t = T (s);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			_M[i] *= t;

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
    operator/           (const S& s) const {

      		assert (cabs(s) != 0.0);
		T t = T (s);
		Matrix<T,P> res = *this;

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] /= t;

		return res;

	}


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator/=         (const Matrix<T,P>& M) {

        size_t i;

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

#if defined USE_VALARRAY
        _M /= M.Container();
#elif defined EXPLICIT_SIMD
        if (fp_type(_M[0]))
        	SSE::binary<T>(_M, M.Container(), SSE::div<T>(), _M);
        else
        	for (size_t i = 0; i < Size(); ++i)
        		_M[i] /= M[i];
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
            _M[i] /= M[i];
#endif

        return *this;

    }


    /**
     * @brief           ELementwise division and assignment operator. i.e. this = this ./ M.
     *
     * @param  M        Divisor matrix.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>&
    operator/=         (const Matrix<S,P> &M) {

        size_t i;

        assert (memcmp(_dim, M.Dim(), dvsz) == 0);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			_M[i] /= M[i];

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
    operator/=         (const S& s) {

    	T zero = T(0.0);
    	T t    = T(s);
        assert (t != zero);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			_M[i] /= T(s);

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
    operator/  (const double s, const Matrix<T,P> &m) {

        Matrix<T,P> res = m;
#ifdef USE_VALARRAY
		res.Container() = s / res.Container();
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < m.Size(); ++i)
            res[i] = s / res[i];
#endif
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
    operator/  (const float s, const Matrix<T,P> &m) {

		Matrix<T,P> res = m;
#ifdef USE_VALARRAY
		res.Container() = s / res.Container();
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < m.Size(); ++i)
            res[i] = s / res[i];
#endif
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
    operator/  (const cxfl s, const Matrix<T,P> &m) {

        Matrix<T,P> res = m;
#ifdef USE_VALARRAY
		res.Container() = s / res.Container();
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < m.Size(); ++i)
            res[i] = s / res[i];
#endif
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
    operator/  (const cxdb s, const Matrix<T,P> &m) {

		Matrix<T,P> res = m;

#ifdef USE_VALARRAY
		res.Container() = s / res.Container();
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < m.Size(); ++i)
            res[i] = s / res[i];
#endif

		return res;

    }


    /**
     * @brief           Elementwise equality with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m == s
     */
    inline friend Matrix<unsigned short>
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
    inline friend Matrix<unsigned short>
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
    inline friend Matrix<unsigned short>
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
    inline friend Matrix<unsigned short>
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
    inline friend Matrix<unsigned short>
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
    inline friend Matrix<unsigned short>
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
    inline friend Matrix<T,P>
    operator&  (const Matrix<unsigned short>& mb, const Matrix<T,P>& m) {
        return   m & mb;
    }

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
        assert (i < INVALID_DIM);
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
        assert (i < INVALID_DIM);
        return _res[i];
    }



    /**
     * @brief           Resolution array
     *
     * @return          All resolutions
     */
    inline const float*
    Res                 () const {
        return &_res[0];
    }



    /**
     * @brief           Get size a given dimension.
     *
     * @param   i       Dimension
     * @return          Dimension
     */
    inline size_t Dim (const size_t i)  const;
    inline const size_t* Dim () const;
    inline const size_t* Dsz () const;
    inline std::vector<size_t> DimVector () const;
    inline size_t Dim (const int i) const;
    inline void Clear ();
    inline const char* GetClassName () const;
    inline void SetClassName (const char* name);
    inline Matrix<T,P> prod (const Matrix<T,P> &M, const char transa, const char transb) const;
    inline Matrix<T,P> prodt (const Matrix<T,P> &M) const;
    inline T dotc (const Matrix<T,P>& M) const;
    inline T dot (const Matrix<T,P>& M) const;
    inline size_t Size () const;
    
};

/*
  

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
*/

#endif // __MATRIX_H__
