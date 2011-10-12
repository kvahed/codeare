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

#include "DWT.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_wavelet2d.h>

Matrix<cplx> 
DWT::Forward (const Matrix<cplx>& m) {

	return Transform (m, false);

}


Matrix<cplx> 
DWT::Backward (const Matrix<cplx>& m) {

	return Transform (m, true);

}


Matrix<cplx> 
DWT::Transform (const Matrix<cplx>& m, const bool bw) {

	Matrix<cplx> res = m;

	gsl_wavelet           *w    = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
	gsl_wavelet_workspace *work = gsl_wavelet_workspace_alloc (res.Size());
	
	double* re = (double*) malloc (res.Size()*sizeof(double));
	double* im = (double*) malloc (res.Size()*sizeof(double));
	size_t  l  = res.Dim(0); // Needs to be power of 2 and quadratic
	
	for (size_t j = 0; j < l; j++)
		for (size_t i = 0; i < l; i++) {
			re [i*l+j] = res.At(j*l+i).real();
			im [i*l+j] = res.At(j*l+i).imag();
		}
	
	if (bw) {

		if (!(gsl_wavelet2d_nstransform_inverse (w, re, l, l, l, work) == GSL_SUCCESS))
			printf ("Wavelet transfor for real part failed\n.");
		if (!(gsl_wavelet2d_nstransform_inverse (w, im, l, l, l, work) == GSL_SUCCESS))
			printf ("Wavelet transfor for imaginary part failed\n.");

	} else {

		if (!(gsl_wavelet2d_nstransform_forward (w, re, l, l, l, work) == GSL_SUCCESS))
			printf ("Wavelet transfor for real part failed\n.");
		if (!(gsl_wavelet2d_nstransform_forward (w, im, l, l, l, work) == GSL_SUCCESS))
			printf ("Wavelet transfor for imaginary part failed\n.");
	}
		
	for (size_t j = 0; j < l; j++)
		for (size_t i = 0; i < l; i++) {
			res.At(j*l+i) = cplx(re[i*l+j],im [i*l+j]);
		}
	
	free (re);
	free (im);

	gsl_wavelet_workspace_free (work);

	return res;

}
