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
   const Vector<size_t>& nx, const NFFT<T>& ft) {

	Matrix< std::complex<T> > out (nx[2],nx[1]);//nk*nc
    for (size_t j = 0; j < nx[1]; ++j)
        Column (out, j, ft * (resize(((nx[0] == 2) ? Slice (sm, j) : Volume (sm, j)),size(in)) * in));
    return out;
}


/**
 * @brief Left hand side operator (i.e. inverse transform) 
 *
 * @param  in      K-space
 * @param  sm      Sensitivities 
 * @param  nx      Sizes & co.
 * @param  ft     FT operators
 * @return         Image
 */
template <class T> inline static Matrix< std::complex<T> >
EH (const Matrix< std::complex<T> >& in, const Matrix< std::complex<T> >& sm,
    const Vector<size_t>& nx, const NFFT<T>& ft) {

	Matrix< std::complex<T> > out = ft ->* Column (in,0) * conj((nx[0] == 2) ? Slice (sm, 0) : Volume (sm, 0));
	for (size_t j = 1; j < nx[1]; ++j)
        out += ft ->* Column (in,j) * conj((nx[0] == 2) ? Slice (sm, j) : Volume (sm, j));
    return out;
}

