/*
 * Vector.hpp
 *
 *  Created on: May 28, 2013
 *      Author: kvahed
 */

#ifndef __VECTOR_HPP__
#define __VECTOR_HPP__

#include "Complex.hpp"

#include <iostream>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>

#ifndef NOEXCEPT
#  ifdef HAVE_CXX11_NOEXCEPT
#    define NOEXCEPT 
//noexcept(true)
#  else
#    define NOEXCEPT
//throw()
#  endif
#endif

#if (_MSC_VER >= 1300) && (WINVER < 0x0500) && !defined(_ftol)
#ifndef _ftol
extern "C" inline long _ftol( double d) { return (long) d;}
#endif
#endif

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
#        include "Allocator.hpp"
#        if defined __AVX__
#            define ALIGNEMENT 32
#        else
#            define ALIGNEMENT 16
#        endif
#        define VECTOR_TYPE(A) std::vector<A,Allocator<A,ALIGNEMENT> >
#        define VECTOR_CONSTR(A,B) std::vector<A,Allocator<A,ALIGNEMENT> >(B)
#        define VECTOR_CONSTR_VAL(A,B,C) std::vector<A,Allocator<A,ALIGNEMENT> >(B,C)
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

template<class T, paradigm P> class Matrix;

/**
 * @brief Alligned data Vector for Matrix<T>
 */
template <class T, paradigm P=SHM> class Vector {
public:

    /**
     * @brief convenience typedefs 
     */
	typedef typename VECTOR_TYPE(T)::iterator iterator;
	typedef typename VECTOR_TYPE(T)::const_iterator const_iterator;

    /**
     * @brief Default constructor
     */
	explicit inline Vector () NOEXCEPT {}

    /**
     * @brief Construct with size
     * @bparam  n  New size
     */
	explicit inline Vector (const size_t n) NOEXCEPT { _data = VECTOR_CONSTR (T,n); }

    /**
     * @brief Construct with size and preset value
     * @param  n  New size
     * @param  val Preset value
     */
	explicit inline Vector (const size_t n, const T& val) NOEXCEPT { _data = VECTOR_CONSTR_VAL(T,n,val); }

    /**
     * @brief Copy constructor from different type
     * @param  cs  To copy
     */
	template<class S> inline Vector (const Vector<S>& cs) NOEXCEPT {
		_data.resize(cs.size());
		for (size_t i = 0; i < _data.size(); ++i)
			_data[i] = (T)cs[i];
	}

#ifdef HAVE_CXX11_RVALUE_REFERENCES
	inline Vector (const Vector<T>& other) NOEXCEPT : _data(other._data) {}
	inline Vector (Vector<T>&& other) NOEXCEPT : _data(std::move(other._data)) {}
	inline Vector& operator= (const Vector<T>& other) NOEXCEPT {
		if (this != &other)
			_data = other._data;
		return *this;
	}
	inline Vector& operator= (Vector<T>&& other) NOEXCEPT {
		if (this != &other)
			_data = std::move(other._data);
		return *this;
	}
#endif

	template<class S> inline Vector<T>& operator= (const Matrix<S,P>& M) {
		*this = M.Container();
		return *this;
	}

	inline Vector<short> operator> (const Vector<T>& v) const {
		assert(v._data.size() == _data.size());
		Vector<short>ret (_data.size());
		for (size_t i = 0; i < _data.size(); ++i)
			ret[i] = CompTraits<T>::greater(_data[i], v[i]);
		return ret;
	}
	inline Vector<short> operator> (const T& t) const {
		Vector<short>ret (_data.size());
		for (size_t i = 0; i < _data.size(); ++i)
			ret[i] = CompTraits<T>::greater(_data[i], t);
		return ret;
	}
	inline Vector<short> operator< (const Vector<T>& v) const {
		assert(v._data.size() == _data.size());
		Vector<short>ret (_data.size());
		for (size_t i = 0; i < _data.size(); ++i)
			ret[i] = CompTraits<T>::less(_data[i], v[i]);
		return ret;
	}
	inline Vector<short> operator< (const T& t) const {
		Vector<short>ret (_data.size());
		for (size_t i = 0; i < _data.size(); ++i)
			ret[i] = CompTraits<T>::less(_data[i], t);
		return ret;
	}
	inline Vector<T> operator& (const Vector<T>& v) const {
		assert(v._data.size() == _data.size());
		Vector<T> ret(_data.size());
		for (size_t i = 0; i < _data.size(); ++i)
			ret[i] = _data[i]&v._data[i];
		return ret;
	}

	/**
     * @brief Elementwise access (lhs)
     * @param  n  n-th element 
     */
	inline T& operator[] (const size_t n) NOEXCEPT { return _data[n]; }
    /**
     * @brief Elementwise access (rhs)
     * @param  n  n-th element 
     */
	inline const T& operator[] (const size_t n) const NOEXCEPT { return _data[n]; }

    /**
     * @brief Access last element (lhs)
     */
	inline T& back () NOEXCEPT { return _data.back(); }
    /**
     * @brief Access last element (rhs)
     */
	inline const T& back () const NOEXCEPT { return _data.back(); }

    /**
     * @brief Access first element (lhs)
     */
	inline T& front () NOEXCEPT { return _data.front(); }
    /**
     * @brief Access first element (rhs)
     */
	inline const T& front () const NOEXCEPT { return _data.front(); }

    /**
     * @brief Access RAM address (lhs)
     * @param  n  at n-th element (default 0)
     */
	inline T* ptr (const size_t n = 0) NOEXCEPT { return &_data[n]; }
    /**
     * @brief Access RAM address (rhs)
     * @param  n  at n-th element (default 0)
     */
	inline const T* ptr (const size_t n = 0) const NOEXCEPT { return &_data[n]; }

    /**
     * @brief Access data vector (lhs)
     */
	inline VECTOR_TYPE(T)& data() NOEXCEPT { return _data; }
    /**
     * @brief Access data vector (rhs)
     */
	inline const VECTOR_TYPE(T)& data() const NOEXCEPT { return _data; }

    /**
     * @brief Vector size
     */
	inline size_t size() const NOEXCEPT { return _data.size(); }

    /**
     * @brief Iterator at start of vector (lhs)
     */
	inline iterator begin() NOEXCEPT { return _data.begin(); }
    /**
     * @brief Iterator at start of vector (rhs)
     */
	inline const_iterator begin() const NOEXCEPT { return _data.begin(); }

    /**
     * @brief Iterator at end of vector (lhs)
     */
	inline iterator end() NOEXCEPT { return _data.end(); }
    /**
     * @brief Iterator at end of vector (rhs)
     */
	inline const_iterator end() const NOEXCEPT { return _data.end(); }

    /**
     * @brief resize data storage
     */
	inline void resize (const size_t n) NOEXCEPT {
		if (!(n==_data.size()))
			_data.resize(n);
	}

    /**
     * @brief resize data storage
     */
	inline void resize (const size_t n, const T val) NOEXCEPT {
		if (!(n==_data.size()))
			_data.resize(n,val);
		else
			_data.assign(n,val);
	}

    /**
     * @brief Add elemet at end
     * @param t  Element to be added
     */
	inline void push_back (const T& t) NOEXCEPT { _data.push_back(t);}
	inline void pop_back () NOEXCEPT { _data.pop_back();}
	inline iterator erase (const iterator& i) NOEXCEPT { return _data.erase(i);}
	inline iterator erase (const iterator& start, const iterator& end) NOEXCEPT { return _data.erase(start, end);}
	inline void insert (const iterator& i, const T& val) NOEXCEPT { _data.insert(i, val);}


	inline void clear() NOEXCEPT {_data.clear();}

	inline bool empty() const NOEXCEPT {return _data.empty();}
	inline bool operator== (const Vector<T>& other) const NOEXCEPT {return _data == other._data;}
	inline bool operator!= (const Vector<T>& other) const NOEXCEPT {return _data != other._data;}
	inline Vector<T>& operator/= (const T& t) NOEXCEPT {
		std::transform (_data.begin(), _data.end(), _data.begin(), std::bind2nd(std::divides<T>(),t));
		return *this;
	}
	inline Vector<T>& operator/= (const Vector<T>& v) NOEXCEPT {
		std::transform (_data.begin(), _data.end(), v.begin(), _data.begin(), std::divides<T>());
		return *this;
	}
	template<class S> inline Vector<T> operator/ (const S& s) const NOEXCEPT {
		Vector<T> ret = *this;
		return ret/=s;
	}
	inline Vector<T>& operator*= (const T& t) NOEXCEPT {
		std::transform (_data.begin(), _data.end(), _data.begin(), std::bind2nd(std::multiplies<T>(),t));
		return *this;
	}
	inline Vector<T>& operator*= (const Vector<T>& v) NOEXCEPT {
		std::transform (_data.begin(), _data.end(), v.begin(), _data.begin(), std::multiplies<T>());
		return *this;
	}
	template<class S> inline Vector<T> operator* (const S& s) const NOEXCEPT {
		Vector<T> ret = *this;
		return ret*=s;
	}
	inline Vector<T>& operator-= (const T& t) NOEXCEPT {
		std::transform (_data.begin(), _data.end(), _data.begin(), std::bind2nd(std::minus<T>(),t));
		return *this;
	}
	inline Vector<T>& operator-= (const Vector<T>& v) NOEXCEPT {
		std::transform (_data.begin(), _data.end(), v.begin(), _data.begin(), std::minus<T>());
		return *this;
	}
	template<class S> inline Vector<T> operator- (const S& s) const NOEXCEPT {
		Vector<T> ret = *this;
		return ret-=s;
	}
	inline Vector<T>& operator+= (const T& t) NOEXCEPT {
		std::transform (_data.begin(), _data.end(), _data.begin(), std::bind2nd(std::plus<T>(),t));
		return *this;
	}
	inline Vector<T>& operator+= (const Vector<T>& v) NOEXCEPT {
		std::transform (_data.begin(), _data.end(), v.begin(), _data.begin(), std::plus<T>());
		return *this;
	}
	template<class S> inline Vector<T> operator+ (const T& s) const NOEXCEPT {
		Vector<T> ret = *this;
		return ret+=s;
	}

private:
	VECTOR_TYPE(T) _data;
};

template<class T> inline static size_t numel (const Vector<T>& v) NOEXCEPT {return v.size();}

template <class T> class vector_inserter {
public:
    Vector<T>& _ct;
    vector_inserter (Vector<T>& ct) NOEXCEPT : _ct(ct) {}
    inline vector_inserter& operator, (const T& val) NOEXCEPT {_ct.push_back(val);return *this;}
};
template <class T> inline vector_inserter<T>& operator+= (Vector<T>& ct,const T& x) NOEXCEPT {
    return vector_inserter<T>(ct),x;
}


template<class T> inline T ct_real (const std::complex<T> ct) NOEXCEPT {return ct.real();}
template<class T> inline T ct_imag (const std::complex<T> ct) NOEXCEPT {return ct.imag();}
template<class T> inline T ct_conj (const T ct) NOEXCEPT {return std::conj(ct);}

template<class T> inline static Vector<T> real (const Vector<std::complex<T> >& c) NOEXCEPT {
	Vector<T> res (c.size());
	std::transform (c.begin(), c.end(), res.begin(), ct_real<T>);
	return res;
}
template<class T> inline static Vector<T> imag (const Vector<std::complex<T> >& c) NOEXCEPT {
	Vector<T> res (c.size());
	std::transform (c.begin(), c.end(), res.begin(), ct_imag<T>);
	return res;
}
template<class T> inline static Vector<T> conj (const Vector<T>& c) NOEXCEPT {
	Vector<T> res (c.size());
	std::transform (c.begin(), c.end(), res.begin(), ct_conj<T>);
	return res;
}
template<class T> inline std::ostream& operator<< (std::ostream& os, const Vector<T>& ct) NOEXCEPT {
	if (ct.size()) {
		for (auto it = ct.begin(); it < ct.end()-1; ++it)
			os << *it << " ";
		os << *(ct.end()-1);
	}
    return os;
}
template<class T> inline static T multiply (const T& a, const T& b) NOEXCEPT {
    return a*b;
}
template<class T> inline static T vprod (const Vector<T>& ct) NOEXCEPT {
	return std::accumulate(ct.begin(), ct.end(), (T)1, multiply<T>);
}
template<class T> inline static T prod (const Vector<T>& ct) NOEXCEPT {
	return std::accumulate(ct.begin(), ct.end(), (T)1, multiply<T>);
}
template<class T> inline static T sum (const Vector<T>& ct) NOEXCEPT {
	return std::accumulate(ct.begin(), ct.end(), (T)0);
}
template<class T> inline static T max (const Vector<T>& ct) NOEXCEPT {
	return *std::max_element(ct.begin(), ct.end());
}
template<class T> inline static T min (const Vector<T>& ct) NOEXCEPT {
	return *std::max_element(ct.begin(), ct.end());
}

template<class T> inline static void swapd (T& x,T& y) NOEXCEPT {T temp=x; x=y; y=temp;}

inline static Vector<size_t> find (const Vector<short>& v) {
	Vector<size_t> ret;
	for (size_t i = 0; i < v.size(); ++i)
		if (v[i])
			ret.push_back(i);
	return ret;
}

#endif /* Vector_HPP_ */
