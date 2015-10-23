/*
 * LocalMaxima.hpp
 *
 *  Created on: Oct 20, 2015
 *      Author: kvahed
 */

#ifndef SRC_MATRIX_CURVES_LOCALMAXIMA_HPP_
#define SRC_MATRIX_CURVES_LOCALMAXIMA_HPP_


#include <Algos.hpp>
#include <Math.hpp>
#include <Print.hpp>

template<class T> struct loc_max_t {
	struct entry {
		size_t pos;
		T val;
		inline entry (size_t p, T v) : pos(p), val(v) {}
	};
	Vector<entry> _db;
	inline void book (const entry& e) {
		_db.PushBack(e);
	}
	inline loc_max_t& operator+= (const entry& e) {
		_db.PushBack(e);
		return *this;
	}
	inline size_t size() const {return _db.size();}
	friend std::ostream& operator<< (std::ostream& os, const loc_max_t& l) {
		std::cout << l.size() << std::endl;
		for (size_t i = 0; i < l.size(); ++i)
			os << l._db[i].pos << " " << l._db[i].val << std::endl;
		return os;
	}
};

//template<class T> inline static loc_max_t<T> findLocalMaxima (const Matrix<T>& A) {
template<class T> inline static Vector<size_t> findLocalMaxima (const Matrix<T>& A) {
	//typedef typename loc_max_t<T>::entry ent;
	//loc_max_t<T> ret;
	Vector<size_t> ret;
	size_t m = size(A,0);
	Matrix<size_t> idx = linspace<size_t>(0,size(A,0)-1,size(A,0));
	Matrix<short> s = sign(diff(A));
	for (size_t i = 1; i < m; ++i)
		if (s[i-1]==1 && s[i]==-1)
			ret.PushBack(i);
	//		ret += ent(i,A[i]);
	return ret;
}
#endif /* SRC_MATRIX_CURVES_LOCALMAXIMA_HPP_ */
