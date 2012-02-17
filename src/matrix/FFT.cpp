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

#include "FFT.hpp"
#include <fftw3.h>


template<> Matrix<cxfl> 
FFT::Forward (const Matrix<cxfl>& m)  {
	
	assert (m.Is1D() || m.Is2D() || m.Is3D());
	
    Matrix<cxfl> res;
	fftwf_plan   p;

	res = Shift(m);

	if (m.Is1D())
		p = fftwf_plan_dft_1d (m.Dim(0),                     (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (m.Is2D())
		p = fftwf_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (m.Is3D())
		p = fftwf_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftwf_execute(p);
	fftwf_destroy_plan(p);

    return Shift(res/sqrt((float)m.Size()));
	
}


template<> Matrix<cxfl>
FFT::Backward (const Matrix<cxfl>& m) {

	assert (m.Is1D() || m.Is2D() || m.Is3D());
	
    Matrix<cxfl> res = Shift(m);
	fftwf_plan   p; 
	
	if (m.Is1D())
		p = fftwf_plan_dft_1d (m.Dim(0),                     (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (m.Is2D())
		p = fftwf_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (m.Is3D())
		p = fftwf_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftwf_complex*)&res[0], (fftwf_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftwf_execute(p);
	fftwf_destroy_plan(p);

	return Shift(res/sqrt((float)m.Size()));
	
}


template<> Matrix<cxdb> 
FFT::Forward (const Matrix<cxdb>& m)  {
	
	assert (m.Is1D() || m.Is2D() || m.Is3D());
	
    Matrix<cxdb> res;
	fftw_plan    p;

	res = Shift(m);

	if (m.Is1D())
		p = fftw_plan_dft_1d (m.Dim(0),                     (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (m.Is2D())
		p = fftw_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	else if (m.Is3D())
		p = fftw_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftw_execute(p);
	fftw_destroy_plan(p);

    return Shift(res / (float)res.Size());
	
}


template<> Matrix<cxdb>
FFT::Backward (const Matrix<cxdb>& m) {

	assert (m.Is1D() || m.Is2D() || m.Is3D());
	
    Matrix<cxdb> res;
	fftw_plan    p; 
	
	res = Shift(m);

	if (m.Is1D())
		p = fftw_plan_dft_1d (m.Dim(0),                     (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (m.Is2D())
		p = fftw_plan_dft_2d (m.Dim(1), m.Dim(0),           (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	else if (m.Is3D())
		p = fftw_plan_dft_3d (m.Dim(2), m.Dim(1), m.Dim(0), (fftw_complex*)&res[0], (fftw_complex*)&res[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_execute(p);
	fftw_destroy_plan(p);

	return Shift(res);
	
}



