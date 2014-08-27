#ifndef __ACCESS_HPP__
#define __ACCESS_HPP__

#include "Algos.hpp"

template<class T> inline size_t
ind2i (const Matrix<T>& M, const size_t& ind) NOEXCEPT {
    return (size_t) ind % size(M,0);
}


template<class T> inline size_t
ind2j (const Matrix<T>& M, const size_t& ind) NOEXCEPT {
    return (size_t) floor (ind/size(M,0)) % (size(M,1)-1);
}


template<class T> inline size_t
ind2k (const Matrix<T>& M, const size_t& ind) NOEXCEPT {
    return (size_t) floor (ind/(size(M,0)*size(M,1))) % (size(M,2)-1);
}


template<class T> inline size_t
ind2l (const Matrix<T>& M, const size_t& ind) NOEXCEPT {
    return (size_t) floor (ind/(size(M,0)*size(M,1)*size(M,2))) % (size(M,3)-1);
}


template<class T> inline size_t
ind2x (const Matrix<T>& M, const size_t& ind, const size_t& dim) NOEXCEPT {

    size_t x = 1;

    for (size_t i = 1; i < dim+1; i++)
        x *= size(M,i-1);

    x = (size_t) floor((double)ind/(double)x) % (size(M,dim));

    return (x);

}


template<class T> inline Matrix<size_t>
ind2sub2d (const Matrix<T>& M, const Matrix<size_t>& inds) NOEXCEPT {

    Matrix<T>      tmp = squeeze(M);
    Matrix<size_t> subs (inds.Size(), 2);

    for(size_t i=0; i < subs.Width(); i++)
        for(size_t j=0; j < subs.Height() ; j++)
            subs(j,i) = ind2x(M,inds(j), i);

    return subs;

}


template <class T> inline Matrix<size_t>
ind2sub3d (const Matrix<T>& M, const Matrix<size_t>& inds) NOEXCEPT {

    Matrix <size_t> subs (inds.Size(), 3);

    for(size_t i=0; i < subs.Width(); i++)
        for(size_t j=0; j < subs.Height() ; j++)
            subs(j,i) = ind2x(M,inds(j), i);

    return subs;

}


template <class T> inline Matrix<size_t>
sub2ind  (const Matrix<T>& M, const Matrix<size_t>& subs) NOEXCEPT {

    size_t n = size(subs,0);

    Matrix<size_t> inds (n);

    /*for (int i = 0; i < n; i++)
      inds[i] = */

    return subs;
}


/**
 * @brief              Get a volume
 * 
 * @param  M           Matrix
 * @param  s           # of volume
 * @return             Desired volume of M
 */
template <class T> inline Matrix<T>
Volume (const Matrix<T>& M, const size_t& s) NOEXCEPT {
	
	Matrix<T> res (size(M,0), size(M,1), size(M,2));
	size_t nc = numel (res);
	
	memcpy (&res[0], M.Ptr(s*nc), nc * sizeof(T));
	
	return res;
	
}
	

/**
 * @brief              Set a volume in a matrix
 * 
 * @param  M           Matrix to insert to
 * @param  s           # of volume
 * @param  A           Matrix to insert
 */
template <class T> inline void
Volume (Matrix<T>& M, const size_t& s, const Matrix<T> A) NOEXCEPT {
	
	assert (size(M,0) == size(A,0));
	assert (size(M,1) == size(A,1));
	assert (size(M,2) == size(A,2));

	size_t nc = size(M,0) * size(M,1) * size(M,2);

	memcpy (&M[s*nc], A.Ptr(), nc*sizeof(T));
	
}
	

/**
 * @brief              Get a slice
 * 
 * @param  M           Matrix
 * @param  s           # of slice
 * @return             Desired slice of M
 */
template <class T> inline Matrix<T>
Slice (const Matrix<T>& M, const size_t& s) NOEXCEPT {
	
	Matrix<T> res (size(M,0),size(M,1));
	size_t nc = numel (res);
	
	memcpy (&res[0], M.Ptr(s*nc), nc*sizeof(T));
	
	return res;
	
}

	
/**
 * @brief              Set a slice
 * 
 * @param  M           Matrix
 * @param  s           # of slice
 * @param  A           New slice
 */
template <class T> inline void
Slice (Matrix<T>& M, const size_t& s, const Matrix<T> A) NOEXCEPT {
	
	assert (size(M,0) == size(A,0));
	assert (size(M,1) == size(A,1));

	size_t ns = size(M,0) * size(M,1);

	memcpy (&M[s * ns], A.Ptr(), ns * sizeof(T));
	
}

	
/**
 * @brief              Set a slice
 * 
 * @param  M           Matrix
 * @param  s           # of slice
 * @param  v           Scalar value           
 */
template <class T> inline void
Slice (Matrix<T>& M, const size_t& s, const T& v) NOEXCEPT {
	size_t ns = size(M,0) * size(M,1);
	Vector<T> vv;
	vv.resize(ns, v);
    std::copy (vv.begin(), vv.end(), M.begin()+s*ns);
}

	
/**
 * @brief              Get a row
 * 
 * @param  M           Matrix
 * @param  r           # of row
 * @return             Desired row of M
 */
template <class T> inline Matrix<T>
Row (const Matrix<T>& M, const size_t& r)  NOEXCEPT {
	
	size_t nc, nr, ns, s, rr;
	Matrix<T> res (1, size(M, 1));

	nc = size(M,1);
	nr = size(M,0);
	ns = nc * nr;

	s  = floor (r/nr);
	rr = r % nr;

	for (size_t i = 0; i < nc; i++)
		res[i] = M[rr + i * nr + s * ns];
	
	return res;
	
}


/**
 * @brief              Copy a row (A) into matrix M
 * 
 * @param  M           Matrix to copy into
 * @param  r           # of row
 * @param  A           Matrix to copy from
 */
template <class T> inline void
Row (Matrix<T>& M, const size_t& r, const Matrix<T> A)  NOEXCEPT {
	
	size_t nc, nr, ns, s, rr;

	nc = size(M,1);
	nr = size(M,0);
	ns = nc * nr;

	assert (size(M,1) == numel (A));
	
	s  = floor (r/nr);
	rr = r % nr;

	for (size_t i = 0; i < nc; i++)
		M[rr + i * nr + s * ns] = A[i];
	
}
	
	
/**
 * @brief              Copy a row (A) into matrix M
 * 
 * @param  M           Matrix to copy into
 * @param  r           # of row
 * @param  v           Scalar value
 */
template <class T> inline void
Row (Matrix<T>& M, const size_t& r, const T& v)  NOEXCEPT {
	
	size_t nc, nr, ns, s, rr;

	nc = size(M,1);
	nr = size(M,0);
	ns = nc * nr;

	s  = floor (r/nr);
	rr = r % nr;

	for (size_t i = 0; i < nc; i++)
		M[rr + i * nr + s * ns] = v;
	
}
	
	
/**
 * @brief              Get a row
 * 
 * @param  M           Matrix
 * @param  c           # of row
 * @return             Desired row of M
 */
template <class T> inline Matrix<T>
Column (const Matrix<T>& M, const size_t& c) NOEXCEPT {
	
	Matrix<T> res (size(M, 0),1);
	
	memcpy (&res[0], M.Ptr(c*size(M, 0)), size(M, 0) * sizeof(T));
	
	return res;
	
}


/**
 * @brief              Get a row
 * 
 * @param  M           Matrix
 * @param  c           # of row
 * @param  A           Vector to copy from
 */
template <class T> inline void
Column (Matrix<T>& M, const size_t& c, const Matrix<T> A) NOEXCEPT {
	
	assert (size(M,0) == size (A,0));
	size_t nc = size(M,0);

	memcpy (&M[c*nc], A.Ptr(), nc * sizeof(T));
	
}


#endif
