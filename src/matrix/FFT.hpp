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
#include "Algos.hpp"

/**
 * @brief 1-3D Discrete Cartesian Fourier transform for Matrix template
 */
	

	/**
	 * @brief    FFT shift
	 *
	 * @param  m To shift
	 * @return   Shifted
	 */
	template <class T> static Matrix<cxfl> 
	fftshift        (const Matrix<T>& m) {

		assert (Is1D(m) || Is2D(m) || Is3D(m));
		
		Matrix<T> res  = m;
		
		for (size_t s = 0; s < m.Dim(2); s++)
			for (size_t l = 0; l < m.Dim(1); l++)
				for (size_t c = 0; c < m.Dim(0); c++)
					res.At (c,l,s) *= (float) pow ((float)-1.0, (float)(s+l+c));
		
		return res;

	}


	template <class T> static Matrix<cxfl> 
	ifftshift        (const Matrix<T>& m) {

		assert (Is1D(m) || Is2D(m) || Is3D(m));
		
		Matrix<T> res  = m;
		
		for (size_t s = 0; s < m.Dim(2); s++)
			for (size_t l = 0; l < m.Dim(1); l++)
				for (size_t c = 0; c < m.Dim(0); c++)
					res.At (c,l,s) *= (float) pow ((float)-1.0, (float)(s+l+c));
		
		return res;

	}

	
	
Matrix<cxfl> 
fft (const Matrix<cxfl>& m)  {
	
	assert (Is1D(m) || Is2D(m) || Is3D(m));
	
    Matrix<cxfl> res;
	fftwf_plan   p;

	res = fftshift(m);

	if      (Is1D(m))
		p = fftwf_plan_dft_1d (m.Dim(0),                     (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (Is2D(m))
		p = fftwf_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (Is3D(m))
		p = fftwf_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftwf_execute(p);
	fftwf_destroy_plan(p);

    return fftshift(res/sqrt((float)m.Size()));
	
}


Matrix<cxfl>
ifft (const Matrix<cxfl>& m) {

	assert (Is1D(m) || Is2D(m) || Is3D(m));
	
    Matrix<cxfl> res = fftshift(m);
	fftwf_plan   p; 
	
	if      (Is1D(m))
		p = fftwf_plan_dft_1d (m.Dim(0),                     (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (Is2D(m))
		p = fftwf_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (Is3D(m))
		p = fftwf_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftwf_execute(p);
	fftwf_destroy_plan(p);

	return fftshift(res/sqrt((float)m.Size()));
	
}


Matrix<cxdb> 
fft (const Matrix<cxdb>& m)  {
	
	assert (Is1D(m) || Is2D(m) || Is3D(m));
	
    Matrix<cxdb> res;
	fftw_plan    p;

	res = fftshift(m);

	if      (Is1D(m))
		p = fftw_plan_dft_1d (m.Dim(0),                     (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (Is2D(m))
		p = fftw_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (Is3D(m))
		p = fftw_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftw_execute(p);
	fftw_destroy_plan(p);

    return fftshift(res / (float)res.Size());
	
}


Matrix<cxdb>
ifft (const Matrix<cxdb>& m) {

	assert (Is1D(m) || Is2D(m) || Is3D(m));
	
    Matrix<cxdb> res;
	fftw_plan    p; 
	
	res = fftshift(m);

	if      (Is1D(m))
		p = fftw_plan_dft_1d (m.Dim(0),                     (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (Is2D(m))
		p = fftw_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (Is3D(m))
		p = fftw_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_execute(p);
	fftw_destroy_plan(p);

	return fftshift(res);
	
}



inline Matrix<double> 
hannwindow (const Matrix<size_t>& size) {

	size_t dim = size.Dim(0);

		assert (dim > 1 && dim < 4);

		Matrix<double> res;

		if      (dim == 1) res = Matrix<double> (size[0], 1);
		else if (dim == 2) res = Matrix<double> (size[0], size[1]);
		else               res = Matrix<double> (size[0], size[1], size[2]);

		float          h, d;
		float          m[3];
		
		if (Is1D(res)) {
			
			m[0] = 0.5 * size[0];
			m[1] = 0.0;
			m[2] = 0.0;
			
		} else if (Is2D(res)) {
			
			m[0] = 0.5 * size[0];
			m[1] = 0.5 * size[1];
			m[2] = 0.0;
			
		} else {
			
			m[0] = 0.5 * size[0];
			m[1] = 0.5 * size[1];
			m[2] = 0.5 * size[2];
			
		}
		
		res = Squeeze(res);
		
		for (size_t s = 0; s < res.Dim(2); s++)
			for (size_t r = 0; r < res.Dim(1); r++)
				for (size_t c = 0; c < res.Dim(0); c++) {
					d = pow( (float)pow(((float)c-m[0])/m[0],2) + pow(((float)r-m[1])/m[1],2) + pow(((float)s-m[2])/m[2],2) , (float)0.5);
					h = (d < 1) ? (0.5 + 0.5 * cos (PI * d)) : 0.0;
					res(c,r,s) = h;
				}
		
		return res;
		
	}


#endif
