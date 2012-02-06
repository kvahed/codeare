/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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

#ifndef __CX_HPP__
#define __CX_HPP__

#include "Matrix.hpp"

class CX {
	
public: 
	
	/**
	 * @brief           Matrix of absolute values  (i.e. |m|)
	 *
	 * @param   m       Input matrix
	 * @return          abs(m)
	 */
	static Matrix<float>
	Abs (const Matrix<cxfl>& m);
	
	/**
	 * @brief           Matrix of absolute values  (i.e. |m|)
	 *
	 * @param   m       Input matrix
	 * @return          abs(m)
	 */
	static Matrix<double>
	Abs (const Matrix<cxdb>& m);
	
	/**
	 * @brief           Matrix of absolute values  (i.e. |m|)
	 *
	 * @param   m       Input matrix
	 * @return          abs(m)
	 */
	static Matrix<float>
	Abs (const Matrix<float>& m);
	
	/**
	 * @brief           Matrix of absolute values  (i.e. |m|)
	 *
	 * @param   m       Input matrix
	 * @return          abs(m)
	 */
	static Matrix<double>
	Abs (const Matrix<double>& m);
	
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<float>
	Arg (const Matrix<cxfl>& m);
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<double>
	Arg (const Matrix<cxdb>& m);
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<float>
	Arg (const Matrix<float>& m);
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<double>
	Arg (const Matrix<double>& m);
	
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<float>
	Real (const Matrix<cxfl>& m);
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<double>
	Real (const Matrix<cxdb>& m);
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<float>
	Real (const Matrix<float>& m);
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<double>
	Real (const Matrix<double>& m);
	
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<float>
	Imag (const Matrix<cxfl>& m);
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<double>
	Imag (const Matrix<cxdb>& m);
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<float>
	Imag (const Matrix<float>& m);
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	static Matrix<double>
	Imag (const Matrix<double>& m);
	
	
	/**
	 * @brief           Matrix of arguments  (i.e. arg(m))
	 *
	 * @param   m       Input matrix
	 * @return          arg(m)
	 */
	template<class T> static Matrix<T>
	Conj (const Matrix<T>& m) {

		Matrix<T> res = m;
		
		if (typeid (T) == typeid (cxfl) || typeid (T) == typeid (cxdb)) {
			
#pragma omp parallel default (shared) 
			{
				
#pragma omp for
				
				for (size_t i = 0; i < m.Size(); i++)
					res[i] = cconj(m[i]);
				
			}		
		
		}

		return res;
	
	}
	
};

#endif

