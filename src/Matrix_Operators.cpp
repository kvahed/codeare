/*template <class T>
Matrix<T> Matrix<T>::operator=(Matrix<T> &M) {
    
    size_t i;
	
    if (nb_alloc) {
		
		if (this->Size() != M.Size()) {
			
			free (_M);
			_M = (T*) malloc (M.Size() * sizeof (T));
			
		}
		
	} else {
		
		_M = (T*) malloc (M.Size() * sizeof (T));
		nb_alloc = 1;
		
	}
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		_dim[i] = M.Dim()[i];
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] = M[i];
		
	}
	
    return *this;

	}*/


template <class T>
Matrix<T> Matrix<T>::operator=(const Matrix<T> &M) {
    
    size_t i;

    if (nb_alloc) {
		
		if (this->Size() != M.Size()) {

			free (_M);
			_M = (T*) malloc (M.Size() * sizeof (T));

		}
		
	} else {

		_M = (T*) malloc (M.Size() * sizeof (T));
		nb_alloc = 1;

	}

	for (size_t i = 0; i < INVALID_DIM; i++)
		_dim[i] = M.Dim()[i];
		
#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())

		for (size_t i = 0; i < Size(); i++)
			_M[i] = M[i];

	}
	
    return *this;

}


template <class T>
Matrix<bool> Matrix<T>::operator==(T s)    {

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
		for (size_t i = 0; i < Size(); i++)
			res[i] = (_M[i] == s);
		
	}
	
    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator>=(T s) {

    Matrix<bool> res(_dim);
    
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] >= s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator<=(T s) {

    Matrix<bool> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] <= s);
    
    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator!=(T s) {

    Matrix<bool> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] != s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator<(T s)    {

    Matrix<bool> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] < s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator>(T s) {

    Matrix<bool> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] > s);

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator->*(Matrix<T> &M) {

    return this->prod(M);

}

template <class T>
Matrix<T> Matrix<T>::operator!() const {

    return this->tr();

}

template <class T>
Matrix<T> Matrix<T>::operator&(Matrix<bool> &M)    {

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

template <class T>
Matrix<T> Matrix<T>::operator&&(Matrix<T> &M) {

    assert(M.Size() == Size());

    Matrix<T> res(_dim);

    for (size_t i = 0; i < Size(); i++)
        res[i] = (M[i] && _M[i]);

    return res;

}

template <class T>
Matrix<T> Matrix<T>::operator||(Matrix<T> M) {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] || M[i]) ? true : false;

	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator==(Matrix<T> M) {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] == M[i]) ? true : false;
	
	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator>=(Matrix<T> M) {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] >= M[i]) ? true : false;

	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator<= (Matrix<T> M) {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] <= M[i]) ? true : false;

	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator!=(Matrix<T> M) {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] != M[i]) ? true : false;

	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator> (Matrix<T> M) {
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] > M[i]) ? true : false;

	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator< (Matrix<T> M) {

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
    for (size_t i = 0; i < Size(); i++)
        res[i] = (_M[i] < M[i]);

	}

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator- () {

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
    for (size_t i = 0; i < Size(); i++)
        res[i] =- _M[i];

	}

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator-(Matrix<T> &M) {

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] - M[i];

	}

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator-(T s) {

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] - s;

	}

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator+(Matrix<T> &M) {

	for (size_t i=0; i < INVALID_DIM; i++)
		assert (Dim(i) == M.Dim(i));

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] + M[i];

	}

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator+(T s) {

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] + s;

	}

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator^(float p) {
    
	Matrix<T> res;
    
	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
	for (size_t i = 0; i < Size(); i++)
        if (p == 0)
            res[i] = 1;
        else
			res[i] = pow(_M[i],  p);


	}

    return res;

}


template <class T> Matrix<T> 
Matrix<T>::operator*(Matrix<T> &M) {

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{

#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		for (size_t i = 0; i < Size(); i++)
			res[i] = _M[i] * M[i];

	}
	
	return res;

}


template <class T> Matrix<T> 
Matrix<T>::operator* (T s) {
    
    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{

#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] * s;
	
	}

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator += (Matrix<T> &M) {
    
    size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] += M[i];
		
	}

    return *this;

}


template <class T> Matrix<T> 
Matrix<T>::operator+= (T s) {
    
    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{

#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] + s;
	
	}

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator *= (Matrix<T> &M) {
    
    size_t i;

	for (size_t i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] *= M[i];
		
	}

    return *this;

}

template <class T>
Matrix<T> Matrix<T>::operator *= (T s) {
    
#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
		for (size_t i = 0; i < Size(); i++)
			_M[i] *= s;
		
	}

    return *this;

}


template <class T> Matrix<T> 
Matrix<T>::operator/(Matrix<T> &M) {

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
	for (size_t i = 0; i < Size(); i++)
		(M[i] != (T)0) ? res[i] = _M[i] / M[i] : 0;

	}

	return res;

}


template <class T> Matrix<T> 
Matrix<T>::operator/ (T s) {
    
	assert (s != (T)0);

    Matrix<T> res;

	for (size_t j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
#pragma omp for schedule (dynamic, Size() / omp_get_num_threads())
		
	for (size_t i = 0; i < Size(); i++)
		res[i] = _M[i] / s;

	}

	return res;

}


template <class T> T           
Matrix<T>::operator[]  (const size_t p) const {
    
    assert(p >= 0);
    assert(p <  Size());
    
    return _M[p];
    
}


template <class T> T           
&Matrix<T>::operator[] (const size_t p) {
    
    assert(p >= 0);
    assert(p <  Size());
    
    return _M[p];
    
}


template <class T>
T Matrix<T>::operator() (const size_t a) const {

    assert(a >= 0);
    assert(a <  Size());

    return _M[a];

}


