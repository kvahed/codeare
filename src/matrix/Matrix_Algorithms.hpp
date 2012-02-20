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

template<class T> inline size_t
Matrix<T>::Ind2i  (const size_t& ind) const { 
	return (size_t) ind % _dim[0];                 
}


template<class T> inline size_t
Matrix<T>::Ind2j  (const size_t& ind) const { 
	return (size_t) floor (ind/_dim[0]) % (_dim[1]-1);
}


template<class T> inline size_t
Matrix<T>::Ind2k  (const size_t& ind) const { 
	return (size_t) floor (ind/(_dim[0]*_dim[1])) % (_dim[2]-1);
}


template<class T> inline size_t
Matrix<T>::Ind2l  (const size_t& ind) const { 
	return (size_t) floor (ind/(_dim[0]*_dim[1]*_dim[2])) % (_dim[3]-1);
}


template<class T> inline size_t
Matrix<T>::Ind2x (const size_t& ind, const size_t& dim) const { 
	
	size_t x = 1;

	for (size_t i = 1; i < dim+1; i++)
		x *= _dim[i-1]; 
	
	x = (size_t) floor((double)ind/(double)x) % (_dim[dim]);

	return (x);

}


template<class T> Matrix<size_t>
Matrix<T>::Ind2Sub2D (const Matrix<size_t>& inds) const {
	
	Matrix<T>      tmp = this->Squeeze();
	Matrix<size_t> subs (inds.Size(), 2);
	
	for(size_t i=0; i < subs.Width(); i++)
		for(size_t j=0; j < subs.Height() ; j++)
			subs(j,i) = Ind2x(inds(j), i);

	return subs; 

}


template <class T> inline Matrix<size_t>
Matrix<T>::Ind2Sub3D (const Matrix<size_t>& inds) const {
	
	Matrix <size_t> subs (inds.Size(), 3);
	
	for(size_t i=0; i < subs.Width(); i++)
		for(size_t j=0; j < subs.Height() ; j++)
			subs(j,i) = Ind2x(inds(j), i);

	return subs; 

}


template <class T> inline Matrix<size_t>
Matrix<T>::Sub2Ind  (const Matrix<size_t>& subs) const {

	size_t n = subs.Dim(0);

	Matrix<size_t> inds (n);

	/*for (int i = 0; i < n; i++)
	  inds[i] = */

	return subs; 
}



