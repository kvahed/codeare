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
#include "DFT.hpp"

/**
 * @brief 1-3D Discrete Cartesian Fourier transform for Matrix template
 */

inline Matrix<cxdb> 
fft (const Matrix<cxdb>& m, const Matrix<double> mask = Matrix<double>(1), const Matrix<cxdb> pc = Matrix<cxdb>(1))  {

	DFT<cxdb> dft (size(m), mask, pc);
    return dft * m;
	
}


inline Matrix<cxdb>
ifft (const Matrix<cxdb>& m, const Matrix<double> mask = Matrix<double>(1), const Matrix<cxdb> pc = Matrix<cxdb>(1)) {
	
	DFT<cxdb> dft (size(m), mask, pc);
    return dft ->* m;
	
}


inline Matrix<cxfl> 
fft (const Matrix<cxfl>& m)  {
	
	DFT<cxfl> dft (size(m));
    return dft * m;
	
}


inline Matrix<cxfl>
ifft (const Matrix<cxfl>& m) {
	
	DFT<cxfl> dft (size(m));
    return dft ->* m;
	
}


#endif
