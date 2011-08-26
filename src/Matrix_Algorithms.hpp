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


template <class T> Matrix<size_t>
Matrix<T>::Ind2Sub2D  (const Matrix<size_t> inds) const {
	
	Matrix<T>      tmp = this->Squeeze();
	Matrix<size_t> subs (inds.Size(), 2);
	
	for(size_t i=0; i < subs.Width(); i++)
		for(size_t j=0; j < subs.Height() ; j++)
			subs(j,i) = Ind2x(inds(j), i);

	return subs; 

}


template <class T> Matrix<size_t>
Matrix<T>::Ind2Sub3D  (const Matrix<size_t> inds) const {
	
	assert(Is2D());

	Matrix <size_t> subs (inds.Size(), 3);
	
	for(size_t i=0; i < subs.Width(); i++)
		for(size_t j=0; j < subs.Height() ; j++)
			subs(j,i) = Ind2x(inds(j), i);

	return subs; 

}


/*
template <class T> Matrix<int>
Sub2Ind  (const Matrix<int>) const {
	matrix < ptrdiff_t > inds(ofMatrixSize(1, subs.width()));
    
	for(size_t i=0; i < subs.width() ; i++) {
		inds(i) = subs(i,subs.height()-1) 
			for(ptrdiff_t j=subs.height()-2 ; j >= 0; j--)
				{
					inds(i) = inds(i)+subs(i, j)*dim(j+1); 
				}
	}
	return subs; 
	} */                                                        



