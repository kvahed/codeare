#include "CX.hpp"

Matrix<double> 
CX::Abs (const Matrix<cxdb>& m) {

	Matrix<double> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
			res[i] = cabs(m[i]);
		
	}		
	
	return res;
	
}
 

Matrix<float> 
CX::Abs (const Matrix<cxfl>& m) {

	Matrix<float> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
			res[i] = cabs(m[i]);
		
	}		
	
	return res;
	
}


Matrix<float> 
CX::Abs (const Matrix<float>& m) {
	
	Matrix<float> res(m);
	return res;
	
}

Matrix<double> 
CX::Abs (const Matrix<double>& m) {
	
	Matrix<double> res(m);
	return res;
	
}


Matrix<double> 
CX::Arg (const Matrix<cxdb>& m) {

	Matrix<double> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
			res[i] = carg(m[i]);
		
	}		
	
	return res;
	
}
 

Matrix<float> 
CX::Arg (const Matrix<cxfl>& m) {

	Matrix<float> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
			res[i] = carg(m[i]);
		
	}		
	
	return res;
	
}


Matrix<float> 
CX::Arg (const Matrix<float>& m) {
	
	Matrix<float> res(m);
	return res;
	
}

Matrix<double> 
CX::Arg (const Matrix<double>& m) {
	
	Matrix<double> res(m);
	return res;
	
}


Matrix<double> 
CX::Real (const Matrix<cxdb>& m) {

	Matrix<double> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
			res[i] = creal(m[i]);
		
	}		
	
	return res;
	
}
 

Matrix<float> 
CX::Real (const Matrix<cxfl>& m) {

	Matrix<float> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
			res[i] = creal(m[i]);
		
	}		
	
	return res;
	
}


Matrix<float> 
CX::Real (const Matrix<float>& m) {
	
	Matrix<float> res(m);
	return res;
	
}

Matrix<double> 
CX::Real (const Matrix<double>& m) {
	
	Matrix<double> res(m);
	return res;
	
}


Matrix<double> 
CX::Imag (const Matrix<cxdb>& m) {

	Matrix<double> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
			res[i] = cimag(m[i]);
		
	}		
	
	return res;
	
}
 

Matrix<float> 
CX::Imag (const Matrix<cxfl>& m) {

	Matrix<float> res(m.Dim());
	
#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		
		for (size_t i = 0; i < m.Size(); i++)
			res[i] = cimag(m[i]);
		
	}		
	
	return res;
	
}


Matrix<float> 
CX::Imag (const Matrix<float>& m) {
	
	Matrix<float> res(m);
	return res;
	
}

Matrix<double> 
CX::Imag (const Matrix<double>& m) {
	
	Matrix<double> res(m);
	return res;
	
}


 


