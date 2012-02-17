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

#ifndef __FFT_HPP__
#define __FFT_HPP__

#include <fftw3.h>
#include "Matrix.hpp"

/**
 * @brief 1-3D Discrete Cartesian Fourier transform for Matrix template
 */
class FFT {
	
public:
	
	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	template <class T> static Matrix<T> 
	Forward      (const Matrix<T>& m);
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	template <class T> static Matrix<T> 
	Backward     (const Matrix<T>& m);
	
	
	/**
	 * @brief    FFT shift
	 *
	 * @param  m To shift
	 * @return   Shifted
	 */
	template <class T> static Matrix<cxfl> 
	Shift        (const Matrix<T>& m) {

		assert (m.Is1D() || m.Is2D() || m.Is3D());
		
		Matrix<T> res  = m;
		
		for (size_t s = 0; s < m.Dim(2); s++)
			for (size_t l = 0; l < m.Dim(1); l++)
				for (size_t c = 0; c < m.Dim(0); c++)
					res.At (c,l,s) *= (float) pow ((float)-1.0, (float)(s+l+c));
		
		return res;

	}
	
	
private:
	
	/**
	 * Static class
	 */
	FFT()  {};
	
	
	/**
	 * Static class
	 */
	~FFT() {};
	
};

#endif
