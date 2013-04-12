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

#ifndef __LINALG_H__
#define __LINALG_H__

template<class T, class S> static int eig (const Matrix<T>& m, Matrix<S>& ev, Matrix<T>& lv, Matrix<T>& rv, const char& jobvl = 'N', const char& jobvr = 'N');
template<class T, class S> static int svd (const Matrix<T>& IN, Matrix<S>& s, Matrix<T>& U, Matrix<T>& V, const char& jobz = 'N');
static Matrix<float> svd (const Matrix<cxfl>& A);
static Matrix<float> svd (const Matrix<float>& A);
static Matrix<double> svd (const Matrix<cxdb>& A);
static Matrix<double> svd (const Matrix<double>& A);
template<class T> static Matrix<T> inv (const Matrix<T>& m); 
template<class T> static Matrix<T> pinv (const Matrix<T>& m, const char& trans = 'N');
template<class T> static Matrix<T> chol (const Matrix<T>& A, const char& uplo = 'U');
template<class T> static double norm (const Matrix<T>& M);
template<class T> static T dotc (const Matrix<T>& A, const Matrix<T>& B);
template<class T> static T DOTC (const Matrix<T>& A, const Matrix<T>& B);
template<class T> T dot (const Matrix<T>& A, const Matrix<T>& B);
template<class T> T DOT  (const Matrix<T>& A, const Matrix<T>& B);
template<class T> static Matrix<T> gemm (const Matrix<T>& A, const Matrix<T>& B, const char& transa = 'N', const char& transb = 'N');
template<class T> Matrix<T> gemv (const Matrix<T>& A, const Matrix<T>& x, const char& trans = 'N');

#endif // __LINALG_H__
