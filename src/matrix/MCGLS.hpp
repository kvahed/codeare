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

class MCGLS {

	
	template<class T> static Matrix<T> 
	Pinv (const Matrix<T>& A, const Matrix<T>& b, const size_t& maxit, const double& conv) {

		// Check x vector 
		// Check A * x valid

		size_t aw = A.Width();

		Matrix<T> AH = !A;
		Matrix<T> x0 = Lapack::GEMV (A, b);
		Matrix<T> x  = Matrix<T>::Zeros (x0.Dim()); 

		return x;

	}

};
