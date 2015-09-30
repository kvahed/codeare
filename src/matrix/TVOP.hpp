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

#ifndef __TVOP_HPP__
#define __TVOP_HPP__

#include "Matrix.hpp"
#include "Operator.hpp"


template<class T, bool b0, bool b1> struct TV2D {};
template<class T> struct TV2D<T,1,1> {
    Matrix<T> inline static fwd (const Matrix<T>& A) {

        size_t M = A.Dim(0), N = A.Dim(1);
        Matrix<T> res (M, N, 2);
                for (size_t n = 0; n < N; ++n)
            for (size_t m = 0; m < M-1; ++m)
                res(m,n,0) = A(m+1,n) - A(m,n);
                for (size_t n = 0; n < N-1; ++n)
            for (size_t m = 0; m < M; ++m)
                res(m,n,1) = A(m,n+1) - A(m,n);
                return res;
    }
    Matrix<T> inline static adj (const Matrix<T>& A) {
        size_t M = A.Dim(0), N = A.Dim(1);
        Matrix<T> res = Matrix<T>(M,N);
        for (size_t m = 0; m < M; ++m)
            for (size_t n = 1; n < N-1; ++n)
                res (n,m) = A(n-1,m,0) - A(n,m,0);
        for (size_t m = 1; m < M-1; ++m)
            for (size_t n = 0; n < N; ++n)
                res (n,m) += A(n,m-1,1) - A(n,m,1);
        return res;
    }
};

template<class T, bool b0, bool b1, bool b2> struct TV3D {};
template<class T> struct TV3D<T,1,1,1> {
    Matrix<T> inline static fwd (const Matrix<T>& A) {
        size_t M = A.Dim(0), N = A.Dim(1), L = A.Dim(2);
        Matrix<T> res = Matrix<T>(M, N, L, 3);
        for (size_t l = 0; l < L; ++l)
            for (size_t n = 0; n < N; ++n)
                for (size_t m = 0; m < M-1; ++m)
                    res(m,n,l,0) = A(m+1,n,l) - A(m,n,l);
        for (size_t l = 0; l < L; ++l)
            for (size_t n = 0; n < N-1; ++n)
                for (size_t m = 0; m < M; ++m)
                    res(m,n,l,1) = A(m,n+1,l) - A(m,n,l);
        for (size_t l = 0; l < L-1; ++l)
            for (size_t n = 0; n < N; ++n)
                for (size_t m = 0; m < M; ++m)
                    res(m,n,l,1) = A(m,n,l+1) - A(m,n,l);
        return res;
    }
    Matrix<T> inline static adj (const Matrix<T>& A) {
        size_t M = A.Dim(0), N = A.Dim(1), L = A.Dim(2);
        Matrix<T> res = Matrix<T>(M,N,L);
        for (size_t l = 0; l < L; ++l)
            for (size_t m = 0; m < M; ++m)
                for (size_t n = 1; n < N-1; ++n)
                    res (n,m,l) = A(n-1,m,l,0) - A(n,m,l,0);
        for (size_t l = 0; l < L; ++l)
            for (size_t m = 1; m < M-1; ++m)
                for (size_t n = 0; n < N; ++n)
                    res (n,m,l) += A(n,m-1,l,1) - A(n,m,l,1);
        for (size_t l = 1; l < L-1; ++l)
            for (size_t m = 0; m < M; ++m)
                for (size_t n = 0; n < N; ++n)
                    res (n,m,l) += A(n,m,l-1,2) - A(n,m,l,2);
        return res;
    }
};

template<class T, bool b0, bool b1, bool b2, bool b3, bool b4> struct TV5D {};
template<class T> struct TV5D<T,0,0,0,1,0> {
    Matrix<T> inline static fwd (const Matrix<T>& A) {
        Vector<size_t> dims = size(A);
        Matrix<T> ret (dims);
        size_t stride = dims[0]*dims[1]*dims[2];
        for (size_t j = 0; j < dims[4]; ++j)
            for (size_t i = 0; i < dims[3]-1; ++i)
                std::transform(A.Begin()+(j*dims[3]+i+1)*stride, A.Begin()+(j*dims[3]+i+2)*stride,
                               A.Begin()+(j*dims[3]+i  )*stride, ret.Begin()+(j*dims[3]+i)*stride,
                               std::minus<T>());
        return ret;
    }
    Matrix<T> inline static adj (const Matrix<T>& A) {
        Vector<size_t> dims = size(A);
        Matrix<T> ret (dims);
        size_t stride = dims[0]*dims[1]*dims[2];
        for (size_t j = 0; j < dims[4]; ++j) {
            size_t offset = j*dims[3];
            for (size_t i = 1; i < dims[3]-1; ++i)
                std::transform(A.Begin()+(offset+i-1)*stride, A.Begin()+(offset+i)*stride,
                               A.Begin()+(offset+i)*stride, ret.Begin()+(offset+i)*stride,
                               std::minus<T>());
        }
        for (size_t j = 0; j < dims[4]; ++j) {
            size_t offset = j*dims[3];
            std::copy_n(A.Begin()+(offset+dims[3]-2)*stride, stride, ret.Begin()+(offset+dims[3]-1)*stride);
            for (size_t i = 0; i < stride; ++i)
                ret[offset*stride+i] = -A[offset*stride+i];
        }
        return ret;
    }
};
template<class T> struct TV5D<T,0,0,0,0,1> {
    Matrix<T> inline static fwd (const Matrix<T>& A) {
        Vector<size_t> dims = size(A);
        Matrix<T> ret (dims);
        size_t stride = dims[0]*dims[1]*dims[2]*dims[3];
        for (size_t i = 0; i < dims[4]-1; i++)
            std::transform(A.Begin()+(i+1)*stride, A.Begin()+(i+2)*stride, A.Begin()+i*stride,
                           ret.Begin()+i*stride, std::minus<T>());
        return ret;
    }
    Matrix<T> inline static adj (const Matrix<T>& A) {
        Vector<size_t> dims = size(A);
        size_t stride = dims[0]*dims[1]*dims[2]*dims[3];
        Matrix<T> ret (dims);
        for (size_t i = 1; i < dims[4]-1; i++)
            std::transform(A.Begin()+(i-1)*stride, A.Begin()+i*stride, A.Begin()+i*stride,
                           ret.Begin()+i*stride, std::minus<T>());
        for (size_t i = 0; i < stride; ++i)
            ret[i] = -A[i];
        std::copy_n(A.Begin()+(dims[4]-2)*stride, stride, ret.Begin()+(dims[4]-1)*stride);
        return ret;
    }
};

enum TVOP_EXCEPTION {UNDEFINED_TV_OPERATOR};

/**
 * @brief 2D Finite difference operator
 */
template <class T>
class TVOP : public Operator<T> {
	

public:

	/**
	 * @brief Default constructor
	 */
	TVOP()  NOEXCEPT {};
    TVOP (const Vector<size_t>& dims) : _dims(dims) {}
    TVOP (const unsigned short dim0, const unsigned short dim1,
          const unsigned short dim2, const unsigned short dim3, const unsigned short dim4) {
        _dims.resize(5);
        _dims[0] = dim0;
        _dims[1] = dim1;
        _dims[2] = dim2;
        _dims[3] = dim3;
        _dims[4] = dim4;
    }


	/**
	 * @brief Default destructor
	 */
	virtual ~TVOP() NOEXCEPT {};


	/**
	 * @brief    2D Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
    
    
    inline Matrix<T> Trafo (const Matrix<T>& A) const {

        Matrix<T> ret;
        
        if      (ndims(A)==2 && (_dims.size()==0 || (_dims[0] == 1 && _dims[1] == 1)))
            ret = TV2D<T,1,1>::fwd(A);
        else if (ndims(A)==3 && (_dims.size()==0 || (_dims[0] == 1 && _dims[1] == 1 && _dims[2] == 1))) 
            ret = TV3D<T,1,1,1>::fwd(A);
        else if (ndims(A)==5 && (_dims.size()==5 && _dims[0] == 0 && _dims[1] == 0 && _dims[2] == 0 && _dims[3] == 1 && _dims[4] == 0))
            ret = TV5D<T,0,0,0,1,0>::fwd(A);
        else if (ndims(A)==5 && (_dims.size()==5 && _dims[0] == 0 && _dims[1] == 0 && _dims[2] == 0 && _dims[3] == 0 && _dims[4] == 1))
            ret = TV5D<T,0,0,0,0,1>::fwd(A);
        else
            throw UNDEFINED_TV_OPERATOR;
        
        return ret;
        
	}	
	
	/**
	 * @brief    Backward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline Matrix<T> Adjoint (const Matrix<T>& A) const {

        Matrix<T> ret;
        
		if      (ndims(A)==3 && (_dims.size()==0 || (_dims[0] == 1 && _dims[1] == 1)))
            ret = TV2D<T,1,1>::adj(A);
		else if (ndims(A)==4 && (_dims.size()==0 || (_dims[0] == 1 && _dims[1] == 1 && _dims[2] == 1)))
            ret = TV3D<T,1,1,1>::adj(A);
		else if (ndims(A)==5 && (_dims.size()==5 && _dims[0] == 0 && _dims[1] == 0 && _dims[2] == 0 && _dims[3] == 1 && _dims[4] == 0))
            ret = TV5D<T,0,0,0,1,0>::adj(A);
        else if (ndims(A)==5 && (_dims.size()==5 && _dims[0] == 0 && _dims[1] == 0 && _dims[2] == 0 && _dims[3] == 0 && _dims[4] == 1))
            ret = TV5D<T,0,0,0,0,1>::adj(A);
        else
            throw UNDEFINED_TV_OPERATOR;
		
		return ret;
		
	}
	


	/**
	 * @brief    Forward transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline Matrix<T> operator* (const Matrix<T>& m) NOEXCEPT {
		return Trafo (m);
	}


	/**
	 * @brief    Adjoint transform
	 *
	 * @param  m To transform
	 * @return   Transform
	 */
	inline Matrix<T> operator->* (const Matrix<T>& m) NOEXCEPT {
		return Adjoint (m);
	}


	inline virtual std::ostream& Print (std::ostream& os) const {
		Operator<T>::Print(os);
		if (_dims.size() == 0)
			os << "    all dims";
		else
			os << "    tv'ed dims: " << _dims;
		return os;
	}

private:
    Vector<unsigned short> _dims;

};

#endif
