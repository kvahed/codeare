#ifndef __ACCESS_HPP__
#define __ACCESS_HPP__

#include "Algos.hpp"

template<class T>
inline size_t ind2i (const Matrix<T>& M, const size_t& ind) NOEXCEPT {
    return (size_t) ind % size(M,0);
}


template<class T>
inline size_t ind2j (const Matrix<T>& M, const size_t& ind) NOEXCEPT {
    return (size_t) floor (ind/size(M,0)) % (size(M,1)-1);
}


template<class T>
inline size_t ind2k (const Matrix<T>& M, const size_t& ind) NOEXCEPT {
    return (size_t) floor (ind/(size(M,0)*size(M,1))) % (size(M,2)-1);
}


template<class T>
inline size_t ind2l (const Matrix<T>& M, const size_t& ind) NOEXCEPT {
    return (size_t) floor (ind/(size(M,0)*size(M,1)*size(M,2))) % (size(M,3)-1);
}


template<class T>
inline size_t ind2x (const Matrix<T>& M, const size_t& ind, const size_t& dim)
		NOEXCEPT {
    size_t x = 1;
    for (size_t i = 1; i < dim+1; i++)
        x *= size(M,i-1);
    x = (size_t) floor((double)ind/(double)x) % (size(M,dim));
    return x;
}

template<class T>
inline Matrix<size_t> ind2sub2d (const Matrix<T>& M,
		const Matrix<size_t>& inds) NOEXCEPT {
    Matrix<T>      tmp = squeeze(M);
    Matrix<size_t> subs (inds.Size(), 2);
    for(size_t i=0; i < subs.Width(); i++)
        for(size_t j=0; j < subs.Height() ; j++)
            subs(j,i) = ind2x(M,inds(j), i);
    return subs;
}

template <class T>
inline Matrix<size_t> ind2sub3d (const Matrix<T>& M, const Matrix<size_t>& inds)
		NOEXCEPT {
    Matrix <size_t> subs (inds.Size(), 3);
    for(size_t i=0; i < subs.Width(); i++)
        for(size_t j=0; j < subs.Height() ; j++)
            subs(j,i) = ind2x(M,inds(j), i);
    return subs;
}

template <class T>
inline Matrix<size_t> sub2ind (const Matrix<T>& M, const Matrix<size_t>& subs)
		NOEXCEPT {
    size_t n = size(subs,0);
    Matrix<size_t> inds (n);
    /*for (int i = 0; i < n; i++)
      inds[i] = */
    return subs;
}

#endif
