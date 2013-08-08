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

#include "OMP.hpp"
#include "Access.hpp"
#include "NFFT.hpp"

/**
 * @brief Right hand side operator E (i.e. forward transform) 
 *
 * @param  in      Image
 * @param  sm      Sensitivities
 * @param  nx      Sizes & co.
 * @param  fts     FT operators
 * @return         K-space 
 */
template <class T> inline static Matrix< std::complex<T> >
E (const Matrix< std::complex<T> >& in, const Matrix< std::complex<T> >& sm,
   const std::vector<size_t>& nx, const std::vector<NFFT<T> >& fts) {

	Matrix< std::complex<T> > out (nx[2],nx[1]);
	
#pragma omp parallel default (shared) 
	{

#pragma omp for 
		for (int j = 0; j < nx[1]; j++)
			Column (out, j, fts[omp_get_thread_num()] * (((nx[0] == 2) ? Slice (sm, j) : Volume (sm, j)) * in));
		
	}
	
	return out;
	
}


/**
 * @brief Left hand side operator (i.e. inverse transform) 
 *
 * @param  in      K-space
 * @param  sm      Sensitivities 
 * @param  nx      Sizes & co.
 * @param  fts     FT operators
 * @return         Image
 */
template <class T> inline static Matrix< std::complex<T> >
EH (const Matrix< std::complex<T> >& in, const Matrix< std::complex<T> >& sm,
    const std::vector<size_t>& nx, const std::vector<NFFT<T> >& fts) {

	Matrix< std::complex<T> > out = zeros< std::complex<T> > (size(sm));

#pragma omp parallel default (shared) 
	{
		
#pragma omp for
		for (int j = 0; j < nx[1]; j++)
			Slice (out, j, fts[omp_get_thread_num()] ->* Column (in,j) * conj(Slice (sm, j)));

	}

	return sum (out, nx[0]);

}

