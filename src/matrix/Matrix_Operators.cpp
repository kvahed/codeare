/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
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

template <class T> inline Matrix<T>&
Matrix<T>::operator= (const Matrix<T>& M) {
    
	if (this->Size() != M.Size())
		_M.resize(M.Size());
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		_dim[i] = M.Dim()[i];
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] = M[i];
		
	}
	
    return *this;
	
}


template <class T> inline Matrix<T> 
Matrix<T>::operator= (const T& s) {
    
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < this->Size(); i++)
			_M[i] = s;
		
	}
	
    return *this;
	
}


template <class T> inline Matrix<bool> 
Matrix<T>::operator== (const T& s) const   {

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			res[i] = (_M[i] == s);
		
	}
	
    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator>= (const T& s) const {

    Matrix<bool> res(_dim);
    
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] >= s);

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator<= (const T& s) const {

    Matrix<bool> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] <= s);
    
    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator!= (const T& s) const {

    Matrix<bool> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] != s);

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator< (const T& s) const {

    Matrix<bool> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] < s);

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator> (const T& s) const {

    Matrix<bool> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] > s);

    return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator->* (const Matrix<T> &M) const {

    return this->prod(M);

}


template <class T> inline Matrix<T> 
Matrix<T>::operator!() const {

    return this->tr();

}


template <class T> inline Matrix<T> 
Matrix<T>::operator& (const Matrix<bool>& M) const {

    for (size_t i = 0; i < INVALID_DIM; i++) 
        assert (_dim[i] == M.Dim()[i]);


    size_t k = 0, i;
    for (i = 0; i < Size(); i++)
        if (M[i])
            k++;
    
    Matrix<T> res(k, 1);
    
    k = 0;
    for (i = 0; i < Size(); i++)
        if (M[i]) {
            res[k] = i;
            k++;
        }

    return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator&& (const Matrix<T>& M) const {

    assert(M.Size() == Size());

    Matrix<T> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (M[i] && _M[i]);

    return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator|| (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] || M[i]) ? true : false;

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator== (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] == M[i]) ? true : false;
	
	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator>= (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] >= M[i]) ? true : false;

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator<= (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] <= M[i]) ? true : false;

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator!= (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] != M[i]) ? true : false;

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator> (const Matrix<T>& M) const {
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] > M[i]) ? true : false;

	}

    return res;

}


template <class T> inline Matrix<bool> 
Matrix<T>::operator< (const Matrix<T>& M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] < M[i]);

	}

    return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator- () const {

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
    for (size_t i = 0; i < Size(); i++)
        res[i] =- _M[i];

	}

    return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator- (const Matrix<S> &M) const {

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] - M[i];

	}

	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator- (const S& s) const {

    Matrix<T> res;
	T t = T(s);

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] - t;

	}

	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator+ (const Matrix<S> &M) const {

	for (size_t i=0; i < INVALID_DIM; i++)
		assert (Dim(i) == M.Dim(i));

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] + M[i];

	}

	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator+ (const S& s) const {

    Matrix<T> res;
	T t = T(s);

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] + t;

	}

	return res;

}


template <class T> inline Matrix<T> 
Matrix<T>::operator^ (const float& p) const {
    
	Matrix<T> res;
    
	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
	for (size_t i = 0; i < Size(); i++)

        if (p == 0)
            res[i] = 1;
        else
			res[i] = pow(_M[i],  p);

	}

    return res;

}


template <class T> inline Matrix<T>
Matrix<T>::operator ^= (const float& p) {
    
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
	for (size_t i = 0; i < Size(); i++)
		
        if (p == 0)
            _M[i] = T(1);
        else
			_M[i] = pow(_M[i],  p);
	
	}

    return *this;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator* (const Matrix<S> &M) const {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<T> res = *this;

#pragma omp parallel default (shared) 
	{

#pragma omp for
		for (size_t i = 0; i < Size(); i++)
			res[i] *= M[i];

	}
	
	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator* (const S& s) const {
    
    Matrix<T> res = *this;
	T t = T(s);

#pragma omp parallel default (shared) 
	{

#pragma omp for
	for (size_t i = 0; i < Size(); i++)
		res[i] *= t;
	
	}

	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator *= (const Matrix<S> &M) {
    
    size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] *= M[i];
		
	}

    return *this;

}


template <class T> template <class S> inline Matrix<T>
Matrix<T>::operator *= (const S& s) {
    
	T t = T(s);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] *= t;
		
	}

    return *this;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator += (const Matrix<S> &M) {
    
    size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] += M[i];
		
	}

    return *this;

}


template <class T> template <class S> inline Matrix<T>
Matrix<T>::operator+= (const S& s) {
    
    Matrix<T> res;
	T t = T(s);

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{

#pragma omp for
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] + t;
	
	}

	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator -= (const Matrix<S>& M) {
    
    size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] -= M[i];
		
	}

    return *this;

}


template <class T> template <class S> inline Matrix<T>
Matrix<T>::operator-= (const S& s) {
    
    Matrix<T> res;
	T t = T(s);

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{

#pragma omp for
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] - t;
	
	}

	return res;

}


template<class T> template<class S> inline Matrix<T> 
Matrix<T>::operator /= (const Matrix<S> &M) {
    
    size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] /= M[i];
		
	}

    return *this;

}


template <class T> template<class S> inline Matrix<T>
Matrix<T>::operator/= (const S& s) {
    
	T t = T(s);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] /= s;
		
	}

    return *this;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator/ (const Matrix<S>& M) const {

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
	for (size_t i = 0; i < Size(); i++)
		(M[i] != (T)0) ? res[i] = _M[i] / M[i] : 0;

	}

	return res;

}


template <class T> template <class S> inline Matrix<T> 
Matrix<T>::operator/ (const S& s) const {
    
	assert (cabs(s) != 0.0);

	T t = T(s);

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for 
		
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] / t;

	}

	return res;

}


template <class T> inline T 
Matrix<T>::operator[]  (const size_t p) const {
    
    assert(p <  Size());
    
    return _M[p];
    
}


template <class T> inline T&
Matrix<T>::operator[] (const size_t p) {
    
    assert(p <  Size());
    
    return _M[p];
    
}


template <class T> inline T 
Matrix<T>::operator() (const size_t a) const {

    assert(a <  Size());

    return _M[a];

}


template<class T> template<class S> inline
Matrix<T>::operator Matrix<S> () const {

	Matrix<S> m (_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (int i = 0; i < this->Size(); i++)
			m[i] = (S)_M[i];
		
	}
	
	return m;

}

