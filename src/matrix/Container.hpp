/*
 * Container.hpp
 *
 *  Created on: May 28, 2013
 *      Author: kvahed
 */

#ifndef CONTAINER_HPP_
#define CONTAINER_HPP_

#include <iostream>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>

#if defined (_MSC_VER) && _MSC_VER<1300
#    define VECTOR_TYPE(A) std::vector<A>
#    define VECTOR_CONSTR(A,B) std::vector<A>(B)
#    define VECTOR_CONSTR_VAL(A,B,C) std::vector<A>(B,C)
#else
#    if defined USE_VALARRAY
#        include <valarray>
#        define VECTOR_TYPE(A) std::valarray<A>
#        define VECTOR_CONSTR(A,B) std::valarray<A>(B)
#    else
#        include "simd/AlignmentAllocator.hpp"
#        if defined __AVX__
#            define ALIGNEMENT 32
#        elif defined __SSE2__
#            define ALIGNEMENT 16
#        endif
#        define VECTOR_TYPE(A) std::vector<A,AlignmentAllocator<A,ALIGNEMENT> >
#        define VECTOR_CONSTR(A,B) std::vector<A,AlignmentAllocator<A,ALIGNEMENT> >(B)
#        define VECTOR_CONSTR_VAL(A,B,C) std::vector<A,AlignmentAllocator<A,ALIGNEMENT> >(B,C)
#    endif
#endif
template<class T> inline static std::ostream&
operator<< (std::ostream& os, const std::vector<T> v) {
    for (size_t i = 0; i < v.size(); ++i)
        os << v[i] << " ";
    return os;
}

typedef unsigned char cbool;

/**
 * @brief   Memory paradigm (share, opencl or message passing)
 */
enum    paradigm {

  SHM, /**< @brief Shared memory (Local RAM) */
  OCL, /**< @brief Open CL GPU RAM */
  MPI  /**< @brief Distributed memory */

};


/**
 * @brief Alligned data container for Matrix<T>
 */
template <class T, paradigm P=SHM> class container {
public:

    /**
     * @brief convenience typedefs 
     */
	typedef typename VECTOR_TYPE(T)::iterator iterator;
	typedef typename VECTOR_TYPE(T)::const_iterator const_iterator;

    /**
     * @brief Default constructor
     */
	explicit inline container () {}
    /**
     * @brief Construct with size
     * @bparam  n  New size
     */
	explicit inline container (const size_t n) { _data = VECTOR_CONSTR (T,n); }
    /**
     * @brief Copy constructor from different type
     * @param  cs  To copy
     */
	template<class S> inline container (const container<S>& cs) {
		_data.resize(cs.size());
		for (size_t i = 0; i < _data.size(); ++i)
			_data[i] = (T)cs[i];
	}

    /**
     * @brief Elementwise access (lhs)
     * @param  n  n-th element 
     */
	inline T& operator[] (const size_t n) { return _data[n]; }
    /**
     * @brief Elementwise access (rhs)
     * @param  n  n-th element 
     */
	inline const T& operator[] (const size_t n) const { return _data[n]; }

    /**
     * @brief Access last element (lhs)
     */
	inline T& back () { return _data.back(); }
    /**
     * @brief Access last element (rhs)
     */
	inline const T& back () const { return _data.back; }

    /**
     * @brief Access first element (lhs)
     */
	inline T& front () { return _data.front(); }
    /**
     * @brief Access first element (rhs)
     */
	inline const T& front () const { return _data.front; }

    /**
     * @brief Access RAM address (lhs)
     * @param  n  at n-th element (default 0)
     */
	inline T* ptr (const size_t n = 0) { return &_data[n]; }
    /**
     * @brief Access RAM address (rhs)
     * @param  n  at n-th element (default 0)
     */
	inline const T* ptr (const size_t n = 0) const { return &_data[n]; }

    /**
     * @brief Access data vector (lhs)
     */
	inline VECTOR_TYPE(T)& data() { return _data; }
    /**
     * @brief Access data vector (rhs)
     */
	inline const VECTOR_TYPE(T)& data() const { return _data; }

    /**
     * @brief Vector size
     */
	inline size_t size() const { return _data.size(); }

    /**
     * @brief Iterator at start of vector (lhs)
     */
	inline iterator begin() { return _data.begin(); }
    /**
     * @brief Iterator at start of vector (rhs)
     */
	inline const_iterator begin() const { return _data.begin(); }

    /**
     * @brief Iterator at end of vector (lhs)
     */
	inline iterator end() { return _data.end(); }
    /**
     * @brief Iterator at end of vector (rhs)
     */
	inline const_iterator end() const { return _data.end(); }

    /**
     * @brief resize data storage
     */
	inline void resize (const size_t n) {
		if (!(n==_data.size()))
			_data.resize(n);
	}

    /**
     * @brief resize data storage
     */
	inline void resize (const size_t n, const T val) {
		if (!(n==_data.size()))
			_data.resize(n,val);
		else
			_data.assign(n,val);
	}

    /**
     * @brief Add elemet at end
     * @param t  Element to be added
     */
	inline void push_back (const T& t) { _data.push_back(t);}

private:
	VECTOR_TYPE(T) _data;
};

template <class T> class vector_inserter {
public:
    container<T>& v;
    vector_inserter(container<T>& v):v(v){}
    vector_inserter& operator,(const T& val){v.push_back(val);return *this;}
};
template <class T> vector_inserter<T>& operator+= (container<T>& v,const T& x){
    return vector_inserter<T>(v),x;
}


template<class T> inline T ct_real (const std::complex<T> ct) {return ct.real();}
template<class T> inline T ct_imag (const std::complex<T> ct) {return ct.imag();}
template<class T> inline T ct_conj (const T ct) {return std::conj(ct);}

template<class T> inline static container<T>
real (const container<std::complex<T> >& c) {
	container<T> res (c.size());
	std::transform (c.begin(), c.end(), res.begin(), ct_real<T>);
	return res;
}
template<class T> inline static container<T>
imag (const container<std::complex<T> >& c) {
	container<T> res (c.size());
	std::transform (c.begin(), c.end(), res.begin(), ct_imag<T>);
	return res;
}
template<class T> inline static container<T>
conj (const container<T>& c) {
	container<T> res (c.size());
	std::transform (c.begin(), c.end(), res.begin(), ct_conj<T>);
	return res;
}
template<class T> inline std::ostream&
operator<< (std::ostream& os, const container<T>& ct) {
    for (typename container<T>::const_iterator it = ct.begin(); it != ct.end(); ++it)
        os << *it << " ";
    return os;
}
template<class T> inline static T multiply (const T& a, const T& b) {
    return a*b;
}
template<class T> inline static T prod (const container<T>& ct) {
	return std::accumulate(ct.begin(), ct.end(), (T)1, multiply<T>);
}
template<class T> inline static T sum (const container<T>& ct) {
	return std::accumulate(ct.begin(), ct.end(), (T)0);
}


#endif /* CONTAINER_HPP_ */
