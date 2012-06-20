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
#include "Algos.hpp"

/**
 * @brief 1-3D Discrete Cartesian Fourier transform for Matrix template
 */
	

inline Matrix<cxfl> 
fft (const Matrix<cxfl>& m)  {
	
	assert (Is1D(m) || Is2D(m) || Is3D(m));
	fftw_import_system_wisdom ();
	
    Matrix<cxfl> res;
	fftwf_plan   p;

	res = fftshift(m);

	if      (Is1D(m))
		p = fftwf_plan_dft_1d (m.Dim(0),                     (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_MEASURE);
	else if (Is2D(m))
		p = fftwf_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_MEASURE);
	else if (Is3D(m))
		p = fftwf_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_MEASURE);
	
	fftwf_execute(p);
	fftwf_destroy_plan(p);

    return fftshift(res/sqrt((float)m.Size()));
	
}


inline Matrix<cxfl>
ifft (const Matrix<cxfl>& m) {

	assert (Is1D(m) || Is2D(m) || Is3D(m));
	fftw_import_system_wisdom ();
	
    Matrix<cxfl> res = fftshift(m);
	fftwf_plan   p; 
	
	if      (Is1D(m))
		p = fftwf_plan_dft_1d (m.Dim(0),                     (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_MEASURE);
	else if (Is2D(m))
		p = fftwf_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_MEASURE);
	else if (Is3D(m))
		p = fftwf_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_MEASURE);
	
	fftwf_execute(p);
	fftwf_destroy_plan(p);

	return fftshift(res/sqrt((float)m.Size()));
	
}


inline Matrix<cxdb> 
fft (const Matrix<cxdb>& m)  {
	
	assert (Is1D(m) || Is2D(m) || Is3D(m));
	
    Matrix<cxdb> res;
	fftw_plan    p;

	res = fftshift(m);

	if      (Is1D(m))
		p = fftw_plan_dft_1d (m.Dim(0),                     (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_MEASURE);
	else if (Is2D(m))
		p = fftw_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_MEASURE);
	else if (Is3D(m))
		p = fftw_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_MEASURE);
	
	fftw_execute(p);
	fftw_destroy_plan(p);

    return fftshift(res / (float)res.Size());
	
}


inline Matrix<cxdb>
ifft (const Matrix<cxdb>& m) {

	assert (Is1D(m) || Is2D(m) || Is3D(m));
	
    Matrix<cxdb> res;
	fftw_plan    p; 
	
	res = fftshift(m);

	if      (Is1D(m))
		p = fftw_plan_dft_1d (m.Dim(0),                     (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_MEASURE);
	else if (Is2D(m))
		p = fftw_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_MEASURE);
	else if (Is3D(m))
		p = fftw_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_MEASURE);
	
	fftw_execute(p);
	fftw_destroy_plan(p);

	return fftshift(res);
	
}





#endif
