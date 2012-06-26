#ifndef __MRI_HPP__
#define __MRI_HPP__


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
		for (size_t i = 0; i < nr; i++) {
			
			for (size_t j = 0; j < nc; j++)
				res[i] += creal(sens(i+j*nr) * conj(sens(i+j*nr)));
			
			res[i] = 1.0 / (((sqroot) ? sqrt (res[i]) : res[i]) + 1.0e-10);
			
		}
		
	}

	return res;

}





#endif
