template <class T>
Matrix<T> Matrix<T>::operator=(Matrix<T> &M) {
    
    int i;

    if (nb_alloc) {
		
		if (this->Size() != M.Size()) {

			delete[](_M);
			_M = (T*) malloc (M.Size() * sizeof (T));

		}
		
	} else {

		_M = (T*) malloc (M.Size() * sizeof (T));
		nb_alloc = 1;

	}

	for (int i = 0; i < INVALID_DIM; i++)
		_dim[i] = M.Dim()[i];
		
#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
		for (int i = 0; i < Size(); i++)
			_M[i] = M[i];
		
	}

    return *this;

}


template <class T>
Matrix<T> Matrix<T>::operator=(const Matrix<T> &M) {
    
    int i;

    if (nb_alloc) {
		
		if (this->Size() != M.Size()) {

			delete[](_M);
			_M = (T*) malloc (M.Size() * sizeof (T));

		}
		
	} else {

		_M = (T*) malloc (M.Size() * sizeof (T));
		nb_alloc = 1;

	}

	for (int i = 0; i < INVALID_DIM; i++)
		_dim[i] = M.Dim()[i];
		
#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)

		for (int i = 0; i < Size(); i++)
			_M[i] = M[i];

	}
	
    return *this;

}


template <class T>
Matrix<bool> Matrix<T>::operator==(T s)    {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] == s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator>=(T s) {

    Matrix<bool> res(_dim);
    
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] >= s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator<=(T s) {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] <= s);
    
    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator!=(T s) {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] != s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator<(T s)    {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] < s);

    return res;

}


template <class T>
Matrix<bool> Matrix<T>::operator>(T s) {

    Matrix<bool> res(_dim);

    for (int i = 0; i < Size(); i++)
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

    for (int i = 0; i < INVALID_DIM; i++) 
        assert (_dim[i] == (int) M.Dim()[i]);


    int k = 0, i;
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

    for (int i = 0; i < Size(); i++)
        res[i] = (M[i] && _M[i]);

    return res;

}

template <class T>
Matrix<T> Matrix<T>::operator||(Matrix<T> M) {

	for (int i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] || M[i]) ? true : false;

	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator==(Matrix<T> M) {

	for (int i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] == M[i]) ? true : false;

	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator>=(Matrix<T> M) {

	for (int i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] >= M[i]) ? true : false;

	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator<= (Matrix<T> M) {

	for (int i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] <= M[i]) ? true : false;

	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator!=(Matrix<T> M) {

	for (int i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] != M[i]) ? true : false;

	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator> (Matrix<T> M) {
	
	for (int i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] > M[i]) ? true : false;

	}

    return res;

}

template <class T>
Matrix<bool> Matrix<T>::operator< (Matrix<T> M) {

	for (int i = 0; i < INVALID_DIM; i++)
		assert (_dim[i] == M.Dim(i));

    Matrix<bool> res(_dim);

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
    for (int i = 0; i < Size(); i++)
        res[i] = (_M[i] < M[i]);

	}

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator- () {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
    for (int i = 0; i < Size(); i++)
        res[i] =- _M[i];

	}

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator-(Matrix<T> &M) {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
	for (int i = 0; i < Size(); i++)
		res[i] = _M[i] - M[i];

	}

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator-(T s) {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
	for (int i = 0; i < Size(); i++)
		res[i] = _M[i] - s;

	}

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator+() {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
    for (int i = 0; i < Size(); i++)
        res[i] =+ _M[i];

	}

    return res;

}


template <class T>
Matrix<T> Matrix<T>::operator+(Matrix<T> &M) {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
	for (int i = 0; i < Size(); i++)
		res[i] = _M[i] + M[i];

	}

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator+(T s) {

    Matrix<T> res;

	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
	for (int i = 0; i < Size(); i++)
		res[i] = _M[i] + s;

	}

	return res;

}


template <class T>
Matrix<T> Matrix<T>::operator^(int p) {
    
	Matrix<T> res;
    
	for (int j = 0; j < INVALID_DIM; j++)
		res.Dim(j) = _dim[j];

	res.Reset();

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = Size() / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
	for (int i = 0; i < Size(); i++)
        if (p == 0)
            res[i] = 1;
        else
            if (p > 0)
                res[i] = power(_M[i], p);
            else
                res[i] = 1 / power(_M[i], -p);

	}

    return res;

}


