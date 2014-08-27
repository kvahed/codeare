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
class TVOP {
	

public:

	/**
	 * @brief    2D Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	template <class T> inline static Matrix<T> 
	Trafo (const Matrix<T>& m) NOEXCEPT {
		
		size_t M = m.Dim(0);
		size_t N = m.Dim(1);
		
		Matrix<T> res (M, N, 2); 
		
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < M-1; ++j)
				res(j,i,0) = m(j+1,i) - m(j,i);
			res (M-1,i,0) = T(0.0);
		}

		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < N-1; ++j)
				res(i,j,1) = m(i,j+1) - m(i,j);
			res (i,N-1,1) = T(0.0);
		}
		
		return res;
		
	}	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	template <class T> inline static Matrix<T> 
	Adjoint    (const Matrix<T>& m) NOEXCEPT {

		size_t M = m.Dim(0);
		size_t N = m.Dim(1);
		
		Matrix<T> resx (M, N);
		Matrix<T> resy (M, N);

		for (int i = 0; i < M; ++i) {
			resy (0,i) = -m(0,i,0);
			for (int j = 1; j < N-1; ++j)
				resy (j,i) = m(j-1,i,0) - m(j,i,0);
			resy (N-1,i) = m(N-2,i,0);
		}

		for (int i = 0; i < N; ++i) {
			resx (i,0) = -m(i,0,1);
			for (int j = 1; j < M-1; ++j)
				resx (i,j) = m(i,j-1,1) - m(i,j,1);
			resx (i,M-1) = m(i,M-2,1);
		}
		
		resx += resy;
		
		return resx;
		
	}
	


	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	template <class T> Matrix<T> 
	operator*    (const Matrix<T>& m) NOEXCEPT {
		return Trafo (m);
	}


	/**
	 * @brief    Adjoint transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	template <class T> Matrix<T> 
	operator->*  (const Matrix<T>& m) NOEXCEPT {
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
