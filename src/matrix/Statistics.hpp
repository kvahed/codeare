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

#ifndef __STATISTICS_HPP__
#define __STATISTICS_HPP__

#include "Lapack.hpp"
#include "Algos.hpp"

/**
 * @brief      Mean reducing a dimension
 *
 * @param  M   Matrix
 * @param  d   Dimension
 * @return     Average of M reducing d matrix
 */
template <class T> static inline Matrix<T>
mean (const Matrix<T>& M, const size_t& d) {
	
	Matrix<T> res  = M;
	float     quot = (float) res.Dim(d);
		
	res = sum (res, d);
	
	return res / quot;
	
}



/**
 * @brief     Covariance of columns
 *
 * @param  m  Matrix of measurements (columns)
 * @return    Covariance
 */
template <class T> static inline Matrix<T> 
cov (const Matrix<T>& m) {

	size_t mm    = size(m,0);
	size_t mn    = size(m,1);

	Matrix<T> mh = mean(m,0);

	Matrix<T> tmp  (mm, mn);

#pragma omp parallel
	{

#pragma omp for schedule (guided)
		for (size_t i = 0; i < mm; i++)
			for (size_t j = 0; j < mn; j++)
				tmp [i*mn+j] = m[i*mn+j] - mh[j];
		
		
	}	

	return (gemm(tmp, tmp, 'C') / (double)(mm-1));

}



#endif // __STATISTICS_HPP__
