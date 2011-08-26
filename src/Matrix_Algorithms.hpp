template<class T> const inline size_t
Matrix<T>::Ind2i  (const size_t& ind) const { 
	return ind % _dim[1];                 
}

template<class T> const inline size_t
Matrix<T>::Ind2j  (const size_t& ind) const { 
	return (ind/_dim[1]) % _dim[2];
}

template<class T> const inline size_t
Matrix<T>::Ind2k  (const size_t& ind) const { 
	return (ind/(_dim[1]*_dim[2])) % _dim[3];
}

template<class T> const inline size_t
Matrix<T>::Ind2l  (const size_t& ind) const { 
	return (ind/(_dim[1]*_dim[2]*_dim[3])) % _dim[4];
}

template<class T> const inline size_t
Matrix<T>::Ind2x  (const size_t& ind) const { 
	
	size_t x = 1;

	for (size_t i = 1; i < ind; i++)
		x *= _dim[i]; 

	return (ind/x) % _dim[ind];

}

/*
template <class T> Matrix<int>
Ind2Sub  (const Matrix<int> inds) const {
	
	Matrix <T> tmp = this->Squeeze();

	Matrix <int> subs (tmp.HDim(), inds.Size());
	
	for(size_t i=0; i < ; i++)
		for(size_t j=0; j <subs.height() ; j++)
			sub(j, i) = ind2subofdim(inds(i), j+1); 
	
	return subs; 

}



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



