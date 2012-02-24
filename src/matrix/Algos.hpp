#ifndef __ALGOS_HPP__
#define __ALGOS_HPP__

template <class T>  inline size_t 
nnz (const Matrix<T>& M) {
	
	size_t nz   = 0;
	T      zero = T(0);
	
	for (int i = 0; i < M.Size(); i++)
		if (M[i] != T(0))
			nz++;
	
	return nz;
	
}


template <class T>  inline bool 
Is1D (const Matrix<T>& M) {
	
	return IsXD(M, 1);
	
}


template <class T>  inline bool 
Is2D (const Matrix<T>& M) {
	
	return IsXD(M, 2);
	
}


template <class T>  inline bool 
Is3D (const Matrix<T>& M) {
	
	return IsXD(M, 3);
	
}


template <class T>  inline bool 
Is4D (const Matrix<T>& M) {
	
	return IsXD(M, 4);
	
}


template <class T>  inline bool 
IsXD (const Matrix<T>& M, const size_t d) {
	
	size_t l = 0;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		if (M.Dim(i) > 1) l++;
	
	return (l == d);
	
}


template <class T>  inline bool 
IsZero (const Matrix<T>& M) {
	
	for (size_t i = 0; i < M.Size(); i++)
		if (M[i] != T(0)) return false;
	
	return true;
	
}


template <class T>  inline bool
IsEmpty (const Matrix<T>& M) {
	
	return (numel(M) == 0);
	
}


template <class T>  inline Matrix<T> 
SOS (Matrix<T>& M, const size_t d) {
	
	assert (M.Dim(d) > 1);
	
	unsigned short nd = HDim(M);
	size_t dim [INVALID_DIM];
		
	for (size_t i = 0; i < INVALID_DIM; i++)
		dim[i] = (i != nd) ? M.Dim(i) : 1;
	
	Matrix<T> res (dim);
		
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < res.Size(); i++) {
			for (size_t j = 0; j < M.Dim(nd); j++)
				res[i] = pow (M[i + j * res.Size()], 2.0);
			pow (res[i], 0.5);
		}
		
	}
	
	return res;
	
}


/*
  template <class T> inline  Matrix<T> 
  SOS (const Matrix<T>& M, const size_t d) {
  
  Matrix<T> res = M;
  
  assert (_dim[d] > 1);
  
  #pragma omp parallel default (shared) 
  {
  
  #pragma omp for
  
  for (int i = 0; i < Size(); i++)
  res[i] = M[i]*M[i];
  
  }
  
  return Sum(res, d);
  
  }
*/

template <class T> inline  Matrix<T>
Squeeze (Matrix<T>& M) {
	
	size_t found = 0;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		if (M.Dim(i) > 1) {
			M.Res(found)   = M.Res(i);
			M.Dim(found++) = M.Dim(i);
		}
	
	for (size_t i = found; i < INVALID_DIM; i++) {
		M.Dim(i) = 1;
		M.Res(i) = 1.0;
	}
	
}


template <class T> inline  Matrix<T>
Mean (const Matrix<T>& M, const size_t d) {
	
	Matrix<T> res  = M;
	float     quot = (float) res.Dim(d);
		
	res = Sum (res, d);
	
	return res / quot;
	
}


template <class T> inline  Matrix<T>
Sum (Matrix<T>& M, const size_t d) {
	
	Matrix<T> res = M;
	
	assert (d < INVALID_DIM);
	
	// No meaningful sum over particular dimension
	if (res.Dim(d) == 1) return res;
	
	// Empty? allocation 
	if (IsEmpty(M))    return res;
	
	// Save old data and resize matrix 
	T* tmp = (T*) malloc (res.Size() * sizeof (T));
	memcpy (tmp, &M[0], res.Size() * sizeof (T));
	
	// Inner size 
	size_t insize = 1;
	for (size_t i = 0; i < d; i++)
		insize *= res.Dim(i);
	
	// Outer size
	size_t outsize = 1;
	for (size_t i = d+1; i < INVALID_DIM; i++)
		outsize *= res.Dim(i);
	
	res.Resize();
	
	// Sum
#pragma omp parallel default (shared) 
	{
		
		size_t tid      = omp_get_thread_num();
		
		for (size_t i = 0; i < outsize; i++) {
			
#pragma omp for
			
			for (size_t j = 0; j < insize; j++) {
				res[j+i*insize] = 0;
				for (size_t k = 0; k < res.Dim(d); k++)
					res[j+i*insize] += tmp[j+i*insize*res.Dim(d)+k*insize];
			}
			
		}
			
	}
	
	// Adjust dminesions and clear tmp
		res.Dim(d) = 1;
		free (tmp);
		
		return res;
		
}


template <class T> inline  size_t
HDim (const Matrix<T>& M) {
	
	size_t nd = 0;
	
	for (size_t i = 0; i < INVALID_DIM; i++)
		nd  = (M.Dim(i) > 1) ? i : nd;
	
	return nd;
	
}


/**
 * @brief           Get the number of matrix cells, i.e. Size * sizeof(T).
 *
 * @param  M        Matrix in question
 * @return          Size in RAM in bytes.
 */
template <class T>  inline size_t
SizeInRAM          (const Matrix<T>& M) {
	
	return Size(M) * sizeof (T);
	
}


/**
 * @brief           Get the number of matrix cells, i.e. dim_0*...*dim_16.
 *
 * @return          Number of cells.
 */
template <class T>  size_t
numel               (const Matrix<T>& M) {
	
	size_t s = 1;
    
	for (size_t i = 0; i < INVALID_DIM; i++)
		s *= M.Dim(i);
    
	return s;
	
}


/**
 * @brief           Get the number of matrix cells, i.e. dim_0*...*dim_16.
 *
 * @return          Number of cells.
 */
template <class T>  Matrix<size_t>
size               (const Matrix<T>& M) {
	
	Matrix<size_t> res (INVALID_DIM, 1);
	size_t ones = 0;
    
	for (size_t i = 0; i < INVALID_DIM; i++) {

		res[i] = M.Dim(i);
		ones   = (res[i] == 1) ? ones + 1 : ones = 0;

	}
    
	return res;
	
}


#endif 
