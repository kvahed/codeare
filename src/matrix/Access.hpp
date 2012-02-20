#ifndef __ACCESS_HPP__
#define __ACCESS_HPP__


#include "Algos.hpp"

class Access {

	
	template <class T> static inline Matrix<T> 
	Volume (const Matrix<T>& M, const size_t s) {
		
		assert (Algos::Is4D(M));
		
		Matrix<T> res;
		
		for (size_t j = 0; j < 3; j++)
			res.Dim(j) = M.Dim(j);
		
		res.Reset();
		
		size_t nc = M.Dim(0)*M.Dim(1)*M.Dim(2);
		
		memcpy (&res[0], &M[s * nc], nc * sizeof(T));
		
		return res;
		
	}
	
	
	template <class T> static inline Matrix<T> 
	Slice (const Matrix<T>& M, const size_t s) {
		
		assert (Algos::Is3D(M));
		
		Matrix<T> res;
		
		for (size_t j = 0; j < 2; j++)
			res.Dim(j) = M.Dim(j);
		
		res.Reset();
		
		size_t nc = M.Dim(0)*M.Dim(1);
		
		memcpy (&res[0], &M[s * nc], nc*sizeof(T));
		
		return res;
		
	}
	
	
	template <class T> static inline Matrix<T> 
	Row (const Matrix<T>& M, const size_t r)  {
		
		assert (Algos::Is2D(M));
		
		Matrix<T> res;
		
		res.Dim(0) = M.Dim(1);
		res.Reset();
		
		for (size_t i = 0; i < M.Dim(1); i++)
			res[i] = M[r + i*M.Dim(0)];
		
		return res;
		
	}
	
	
	template <class T> static inline Matrix<T> 
	Column (const Matrix<T>& M, const size_t c) {
		
		Matrix<T> res;
		
		res.Dim(0) = M.Dim(0);
		res.Reset();
		
		memcpy (&res[0], M[c*M.Dim(0)], M.Dim(0) * sizeof(T));
		
		return res;
		
	}
	
};
