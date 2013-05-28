/*
 * Container.hpp
 *
 *  Created on: May 28, 2013
 *      Author: kvahed
 */

#ifndef CONTAINER_HPP_
#define CONTAINER_HPP_

#include "SIMD.hpp"

#include <assert.h>

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
	inline container (const size_t n) { _data = VECTOR_CONSTR (T,n); }
	inline container () { _data = VECTOR_CONSTR (T,1); }
	inline T& operator[] (const size_t n) { return _data[n]; }
	inline T operator[] (const size_t n) const { return _data[n]; }
	inline const T* memory (const size_t n) const { return &_data[n]; }
	inline T* memory (const size_t n)  { return &_data[n]; }
	inline size_t size() const { return _data.size(); }
	inline ~container () {}
	typedef typename VECTOR_TYPE(T)::iterator iterator;
	typedef typename VECTOR_TYPE(T)::const_iterator const_iterator;
	inline iterator begin() { return _data.begin(); }
	inline iterator end() { return _data.end(); }
	inline const_iterator begin() const { return _data.begin(); }
	inline const_iterator end() const { return _data.end(); }
	inline void resize (const size_t n, const T val = T()) { _data.resize(n,val); }
private:
	VECTOR_TYPE(T) _data;
};

template<> template<>
class container<bool,SHM> {
public:
	inline container (const size_t n) { _n = n; _data = (bool*) malloc (n * sizeof(bool)); _malloc = true;}
	inline container () { _n = 1; _data = (bool*) malloc (sizeof(bool)); _malloc = true;}
	inline bool& operator[] (const size_t n) { assert (n<_n); return _data[n]; }
	inline bool operator[] (const size_t n) const { assert (n<_n); return _data[n]; }
	inline const bool* memory (const size_t n) const { return &_data[n]; }
	inline bool* memory (const size_t n)  { return &_data[n]; }
	inline size_t size() const { return _n; }
	inline ~container () {	/*if (_malloc) delete[] _data;*/ }
	inline void resize (const size_t n, const bool val = bool()) {
		_data = (bool*) realloc (_data, n * sizeof(bool)); _malloc = true;}
private:
	bool *_data;
	size_t _n;
	bool _malloc;
};


#endif /* CONTAINER_HPP_ */
