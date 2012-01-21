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
	static Matrix<cplx> 
	Forward      (const Matrix<cplx>& m);
	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline static Matrix<cplx> 
	Backward     (const Matrix<cplx>& m) {

		assert (m.Is1D() || m.Is2D() || m.Is3D());
		
		Matrix<cplx> res = m;
		fftwf_plan   p; 
	
		res = Shift(res);
		
		if (m.Is1D())
			p = fftwf_plan_dft_1d (m.Dim(0),                     (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
		else if (m.Is2D())
			p = fftwf_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
		else if (m.Is3D())
			p = fftwf_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
		
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		
		res = Shift(res);
		
		return res;
		
	};
	
	
	/**
	 * @brief    FFT shift
	 *
	 * @param  m To shift
	 * @return   Shifted
	 */
	static Matrix<cplx> 
	Shift        (const Matrix<cplx>& m);
	
	
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
