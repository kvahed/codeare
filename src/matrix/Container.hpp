/*
 * Container.hpp
 *
 *  Created on: May 28, 2013
 *      Author: kvahed
 */

#ifndef CONTAINER_HPP_
#define CONTAINER_HPP_

#include "SIMD.hpp"

#include <iostream>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <boost/dynamic_bitset.hpp>

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
	inline container () { _data = VECTOR_CONSTR (T,1); }
	inline container (const size_t n) { _data = VECTOR_CONSTR (T,n); }
	inline T& operator[] (const size_t n) { return _data[n]; }
	inline T operator[] (const size_t n) const { return _data[n]; }
	inline const T* memory (const size_t n = 0) const { return &_data[n]; }
	inline T* memory (const size_t n = 0)  { return &_data[n]; }
	inline VECTOR_TYPE(T) data() const { return _data; }
	inline VECTOR_TYPE(T)& data() { return _data; }
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
/*
template<> template<>
class container<bool> {
public:
	inline container () { _data = boost::dynamic_bitset<>(1);}
	inline container (const size_t n) { _data = boost::dynamic_bitset<>	(n); }
	inline boost::dynamic_bitset<>::reference& operator[] (const size_t n) { return _data[n]; }
	inline bool operator[] (const size_t n) const { return _data[n]; }
	inline boost::dynamic_bitset<> data() const { return _data; }
	inline boost::dynamic_bitset<>& data() { return _data; }
	inline size_t size() const { return _data.size(); }
	inline ~container () {}
	inline void resize (const size_t n, const bool val = bool()) { _data.resize(n,val); }
private:
	boost::dynamic_bitset<> _data;
};
*/

template<> template<>
class container<bool,SHM> {
public:
	inline container (const size_t n) : _malloc(false), _data(NULL), _n(0) {
		_n = n; _data = (bool*) malloc (ceil(_n/8)); _malloc = true;
	}
	inline container () : _malloc(false), _data(NULL), _n(0) {
		_n = 1; _data = (bool*) malloc (1); _malloc = true;
	}
	inline bool& operator[] (const size_t n) { return _data[n]; }
	inline bool operator[] (const size_t n) const { assert (n<_n); return _data[n]; }
	inline const bool* memory (const size_t n = 0) const { return &_data[n]; }
	inline bool* memory (const size_t n = 0)  { return &_data[n]; }
	inline const bool* data() const { return _data; }
	inline bool* data() { return _data; }
	inline size_t size() const { return _n; }
	inline ~container () { /*if (_malloc) _data = (bool*) realloc(_data,0); _malloc = false;*/ }
	inline void resize (const size_t n, const bool val = bool()) {
		_n = n;
		_data = (_malloc) ?
			(bool*) realloc(_data, ceil(_n/8)) : (bool*) malloc (ceil(_n/8));
		_malloc = true;
	}
	inline container<bool,SHM> operator= (const container<bool,SHM>& ct) {
		_n = ct.size();
		if (this != &ct) {
			_data = (_malloc) ?
				(bool*) realloc(_data, ceil(_n/8)) : (bool*) malloc (ceil(_n/8));
			for (size_t i = 0; i < _n; i++)
				_data[i] = ct[i];
		}
		return *this;
	}
private:
	bool *_data;
	size_t _n;
	bool _malloc;
};


#endif /* CONTAINER_HPP_ */
