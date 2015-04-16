/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301  USA
 */

#ifndef __TVOP_HPP__
#define __TVOP_HPP__

#include "Matrix.hpp"


/**
 * @brief 2D Finite difference operator
 */
template <class T>
class TVOP {
	

public:

	/**
	 * @brief    2D Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline static Matrix<T> Trafo (const Matrix<T>& A) NOEXCEPT {
		
		Matrix<T> res;

		if (ndims(A)==2) {
			size_t M = A.Dim(0), N = A.Dim(1);
			res = Matrix<T>(M, N, 2);
			for (size_t n = 0; n < N; ++n)
				for (size_t m = 0; m < M-1; ++m)
					res(m,n,0) = A(m+1,n) - A(m,n);
			for (size_t n = 0; n < N-1; ++n)
				for (size_t m = 0; m < M; ++m)
					res(m,n,1) = A(m,n+1) - A(m,n);
		} else if (ndims(A)==3) {
			size_t M = A.Dim(0), N = A.Dim(1), L = A.Dim(2);
			res = Matrix<T>(M, N, L, 3);
			for (size_t l = 0; l < L; ++l)
				for (size_t n = 0; n < N; ++n)
					for (size_t m = 0; m < M-1; ++m)
						res(m,n,l,0) = A(m+1,n,l) - A(m,n,l);
			for (size_t l = 0; l < L; ++l)
				for (size_t n = 0; n < N-1; ++n)
					for (size_t m = 0; m < M; ++m)
						res(m,n,l,1) = A(m,n+1,l) - A(m,n,l);
			for (size_t l = 0; l < L-1; ++l)
				for (size_t n = 0; n < N; ++n)
					for (size_t m = 0; m < M; ++m)
						res(m,n,l,1) = A(m,n,l+1) - A(m,n,l);
		}
		
		return res;
		
	}	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline static Matrix<T> Adjoint (const Matrix<T>& A) NOEXCEPT {

		Matrix<T> resx, resy;

		if (ndims(A)==3) {
			size_t M = A.Dim(0), N = A.Dim(1);
			resx = Matrix<T>(M,N);
			resy = Matrix<T>(M,N);
			for (size_t m = 0; m < M; ++m)
				for (size_t n = 1; n < N-1; ++n)
					resy (n,m) = A(n-1,m,0) - A(n,m,0);
			for (size_t m = 1; m < M-1; ++m)
				for (size_t n = 0; n < N; ++n)
					resx (n,m) = A(n,m-1,1) - A(n,m,1);
			resx += resy;
		} else if (ndims(A)==4) {
			size_t M = A.Dim(0), N = A.Dim(1), L = A.Dim(2);
			resx = Matrix<T>(M,N,L);
			resy = Matrix<T>(M,N,L);
			Matrix<T> resz = Matrix<T>(M,N,L);
			for (size_t l = 0; l < L; ++l)
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 1; n < N-1; ++n)
						resy (n,m,l) = A(n-1,m,l,0) - A(n,m,l,0);
			for (size_t l = 0; l < L; ++l)
				for (size_t m = 1; m < M-1; ++m)
					for (size_t n = 0; n < N; ++n)
						resx (n,m,l) = A(n,m-1,l,1) - A(n,m,l,1);
			for (size_t l = 1; l < L-1; ++l)
				for (size_t m = 0; m < M; ++m)
					for (size_t n = 0; n < N; ++n)
						resz (n,m,l) = A(n,m,l-1,2) - A(n,m,l,2);
			resx += resy + resz;
		}
		
		return resx;
		
	}
	


	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline Matrix<T> operator* (const Matrix<T>& m) NOEXCEPT {
		return Trafo (m);
	}


	/**
	 * @brief    Adjoint transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline Matrix<T> operator->* (const Matrix<T>& m) NOEXCEPT {
		return Adjoint (m);
	}


	/**
	 * @brief Default constructor
	 */
	TVOP()  NOEXCEPT {};


	/**
	 * @brief Default destructor
	 */
	~TVOP() NOEXCEPT {};

};

#endif
