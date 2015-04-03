/*
 *  codeare Copyright (C) 2007-2015 Kaveh Vahedipour
 *                                  Forschungszentrum Juelich, Germany
 *                                  NYU School of Medicine, New York, USA
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
   const Vector<size_t>& nx, const Vector<NFFT<T> >& ft) {
	Matrix< std::complex<T> > out (nx[2],nx[1]);
#pragma omp parallel for schedule (guided, 1)
    for (int j = 0; j < nx[1]; ++j) {
        int k = omp_get_thread_num();
        Column (out, j, ft[k] * (resize(((nx[0] == 2) ? Slice (sm, j) : Volume (sm, j)),size(in)) * in));
    }
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
    const Vector<size_t>& nx, const Vector<NFFT<T> >& ft) {
    Vector<size_t> n = size((nx[0] == 2) ? Slice (sm, 0) : Volume (sm, 0));
    n.PushBack(nx[1]);
	Matrix< std::complex<T> > out (n);
#pragma omp parallel for schedule (guided, 1)
	for (int j = 0; j < nx[1]; ++j) {
        int k = omp_get_thread_num();
        if (nx[0] == 2)
            Slice  (out, j, ft[k] ->* Column (in,j) * conj(Slice  (sm, j)));
        else
            Volume (out, j, ft[k] ->* Column (in,j) * conj(Volume (sm, j)));
    }
    return squeeze(sum(out,n.size()-1));
}
