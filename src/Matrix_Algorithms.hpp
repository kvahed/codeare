template<class T> const inline size_t
Matrix<T>::Ind2i  (const size_t& ind) const { 
	return (size_t) ind % _dim[0];                 
}



template<class T> const inline size_t
Matrix<T>::Ind2j  (const size_t& ind) const { 
	return (size_t) floor (ind/_dim[0]) % (_dim[1]-1);
}



template<class T> const inline size_t
Matrix<T>::Ind2k  (const size_t& ind) const { 
	return (size_t) floor (ind/(_dim[0]*_dim[1])) % (_dim[2]-1);
}



template<class T> const inline size_t
Matrix<T>::Ind2l  (const size_t& ind) const { 
	return (size_t) floor (ind/(_dim[0]*_dim[1]*_dim[2])) % (_dim[3]-1);
}


template<class T> inline const size_t
Matrix<T>::Ind2x (const size_t& ind, const size_t& dim) const { 
	
	size_t x = 1;

	for (size_t i = 1; i < dim+1; i++)
		x *= _dim[i-1]; 
	
	x = (size_t) floor((double)ind/(double)x) % (_dim[dim]);

	return (x);

}


template<class T> inline Matrix<size_t>
Matrix<T>::Ind2Sub2D (const Matrix<size_t> inds) const {
	
	Matrix<T>      tmp = this->Squeeze();
	Matrix<size_t> subs (inds.Size(), 2);
	
	for(size_t i=0; i < subs.Width(); i++)
		for(size_t j=0; j < subs.Height() ; j++)
			subs(j,i) = Ind2x(inds(j), i);

	return subs; 

}


template <class T> inline Matrix<size_t>
Matrix<T>::Ind2Sub3D (const Matrix<size_t> inds) const {
	
	assert(Is2D());

	Matrix <size_t> subs (inds.Size(), 3);
	
	for(size_t i=0; i < subs.Width(); i++)
		for(size_t j=0; j < subs.Height() ; j++)
			subs(j,i) = Ind2x(inds(j), i);

	return subs; 

}


template <class T> inline Matrix<size_t>
Matrix<T>::Sub2Ind  (const Matrix<size_t> subs) const {

	int n = subs.Dim(0);

	Matrix<size_t> inds (n);

	/*for (int i = 0; i < n; i++)
	  inds[i] = */

	return subs; 
}


template<class T> inline const bool 
Matrix<T>::Empty () const {

	_M.size() == 0;

}


template <class T> inline const bool 
Matrix<T>::IsXD (const size_t d) const {

	size_t l = 0;

	for (size_t i = 0; i < INVALID_DIM; i++)
		if (_dim[i] > 1) l++;

	return (l == d);

}


template <class T> inline const bool 
Matrix<T>::Is1D () const {
	
	return IsXD(1);

}


template <class T> inline const bool 
Matrix<T>::Is2D () const {
	
	return IsXD(2);

}


template <class T> inline const bool 
Matrix<T>::Is3D () const {
	
	return IsXD(3);

}


template <class T> inline const bool 
Matrix<T>::Is4D () const {
	
	return IsXD(4);

}



template<class T> Matrix<T>
Matrix<T>::Resample (float f, Interpolation i) {

	// We handle only isotropic resampling.
	assert ((float)(_dim[0] / (size_t)f) == (float)_dim[0] / f);

	Matrix<T> res(this->Dim(0)/f, this->Dim(1)/f, this->Dim(2)/f);
	

	return res;

}
