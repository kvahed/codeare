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
#        include "AlignmentAllocator.hpp"
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



template <class T, paradigm P=SHM>
class container {
public:
	typedef typename VECTOR_TYPE(T)::iterator iterator;
	typedef typename VECTOR_TYPE(T)::const_iterator const_iterator;
	inline container () { _data = VECTOR_CONSTR (T,1); }
	inline container (const size_t n) { assert(n>0); _data = VECTOR_CONSTR (T,n); }
	inline T& operator[] (const size_t n) { return _data[n]; }
	inline T operator[] (const size_t n) const { return _data[n]; }
	inline T& back () { return _data.back(); }
	inline T back () const { return _data.back; }
	inline T& front () { return _data.front(); }
	inline T front () const { return _data.front; }
	inline T& at (const size_t n) { return _data.at(n); }
	inline T at (const size_t n) const { return _data.at(n); }
	inline const T* memory (const size_t n = 0) const { return &_data[n]; }
	inline T* memory (const size_t n = 0)  { return &_data[n]; }
	inline VECTOR_TYPE(T) data() const { return _data; }
	inline VECTOR_TYPE(T)& data() { return _data; }
	inline size_t size() const { return _data.size(); }
	inline ~container () {}
	inline iterator begin() { return _data.begin(); }
	inline iterator end() { return _data.end(); }
	inline const_iterator begin() const { return _data.begin(); }
	inline const_iterator end() const { return _data.end(); }
	inline void resize (const size_t n, const T val = T()) { assert(n>0); _data.resize(n,val); }
private:
	VECTOR_TYPE(T) _data;
};

template<class T>
inline T ct_real (const std::complex<T> ct) {return ct.real();};
template<class T>
inline T ct_imag (const std::complex<T> ct) {return ct.imag();};
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

template<class T> inline std::ostream&
operator<< (std::ostream& os, const container<T>& ct) {
    for (typename container<T>::const_iterator it = ct.begin(); it != ct.end(); ++it)
        os << *it << " ";
    return os;
}

#endif /* CONTAINER_HPP_ */
