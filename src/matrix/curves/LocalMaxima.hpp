/*
 * LocalMaxima.hpp
 *
 *  Created on: Oct 20, 2015
 *      Author: kvahed
 */

#ifndef SRC_MATRIX_CURVES_LOCALMAXIMA_HPP_
#define SRC_MATRIX_CURVES_LOCALMAXIMA_HPP_

/*
function [iPk, iInflect] = findLocalMaxima(yTemp)
% bookend Y by NaN and make index vector
yTemp = [NaN; yTemp; NaN];
iTemp = (1:numel(yTemp)).';

% keep only the first of any adjacent pairs of equal values (including NaN).
yFinite = ~isnan(yTemp);
iNeq = [1; 1 + find((yTemp(1:end-1) ~= yTemp(2:end)) & ...
                    (yFinite(1:end-1) | yFinite(2:end)))];
iTemp = iTemp(iNeq);

% take the sign of the first sample derivative
s = sign(diff(yTemp(iTemp)));

% find local maxima
iMax = 1 + find(diff(s)<0);

% find all transitions from rising to falling or to NaN
iAny = 1 + find(s(1:end-1)~=s(2:end));

% index into the original index vector without the NaN bookend.
iInflect = iTemp(iAny)-1;
iPk = iTemp(iMax)-1;
*/
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

template<class T> inline static loc_max_t<T> findLocalMaxima (const Matrix<T>& A) {
	typedef typename loc_max_t<T>::entry ent;
	loc_max_t<T> ret;
	size_t m = size(A,0);
	Matrix<size_t> idx = linspace<size_t>(0,size(A,0)-1,size(A,0));
	Matrix<short> s = sign(diff(A));
	for (size_t i = 1; i < m; ++i)
		if (s[i-1]==1 && s[i]==-1)
			ret += ent(i,A[i]);
	return ret;
}
#endif /* SRC_MATRIX_CURVES_LOCALMAXIMA_HPP_ */
