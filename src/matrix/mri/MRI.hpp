#ifndef __MRI_HPP__
#define __MRI_HPP__

#include "Algos.hpp"

template <class T> inline static Matrix<T>
IntensityMap (const Matrix< std::complex <T> >& sens, bool sqroot = true) {

	size_t dim = ndims(sens)-1;
	size_t nc  = size(sens,dim);
	size_t nr  = numel(sens)/nc;
	
	Matrix<size_t> dims = size (sens);
	dims [dim] = 1;

	Matrix<T> res = zeros<T> (dims);
	
#pragma omp parallel default (shared)
	{		
#pragma omp for schedule (guided)
		for (int i = 0; i < nr; i++) {
			for (size_t j = 0; j < nc; j++)
				res[i] += real(sens(i+j*nr) * TypeTraits<std::complex<T> >::Conj(sens(i+j*nr)));
			res[i] = (sqroot) ?
				1. / (sqrt(res[i]) + 1.e-9):
				1. / (     res[i]  + 1.e-9);
		}
	}
	return squeeze(res);
}

template <class T> inline static Matrix<T>
phase_combine (const Matrix<T>& M, const size_t d) {
    
	Matrix<size_t> sz = size(M);
	assert (d < sz.Size());
	size_t        dim = sz[d];
	assert (dim == 2);
	Matrix<T>     ret;

	// Empty? allocation
	if (isempty(M))
		return ret;

	// Inner size
	size_t insize = 1;
	for (size_t i = 0; i < d; ++i)
		insize *= sz[i];

	// Outer size
	size_t outsize = 1;
	for (size_t i = d+1; i < MIN(M.NDim(),numel(sz)); ++i)
		outsize *= sz[i];

	// Adjust size vector and allocate
	sz [d] = 1;
	ret = zeros<T>(sz);

	// Combine
#pragma omp parallel default (shared)
	{
#pragma omp for

		for (size_t i = 0; i < outsize; ++i)
			for (size_t j = 0; j < insize; ++j)
				ret[i*insize + j] += sqrt(M[j] * TypeTraits<T>::Conj(M[i*insize + j]));

	}

	return ret;
}


#endif
