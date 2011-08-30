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

#include "Lapack.hpp"

template <class T>
int  Matrix<T>::EIG (Matrix<T>* ev, Matrix<T>* lev, Matrix<T>* rev, const bool cv) {

	// 2D square matrix.
	if (!Is2D())
		return -2;
	if (_dim[COL] != _dim[LIN])
		return -3;
	
	char   jobvl = (cv) ? 'V' : 'N';
	char   jobvr = (cv) ? 'V' : 'N';
	
	int    n     = _dim[COL];
	int    lda   = n;
    int    ldvl  =  (cv) ? n : 1;
    int    ldvr  =  (cv) ? n : 1;
	int    info  =  0;
	int    lwork = -1;
	
	T*     wi    = (T*)     malloc (  n * sizeof(T));
	T*     work  = (T*)     malloc (  1 * sizeof(T));
	float* rwork = (float*) malloc (2*n * sizeof(float));

	// Workspace query
	if (typeid(T) == typeid(cplx))
		cgeev_ (&jobvl, &jobvr, &n, &this->At(0), &lda, &ev->At(0),     &lev->At(0), &ldvl, &rev->At(0), &ldvr, work, &lwork, rwork, &info);
	else if (typeid(T) == typeid(double))
		dgeev_ (&jobvl, &jobvr, &n, &this->At(0), &lda, &ev->At(0), wi, &lev->At(0), &ldvl, &rev->At(0), &ldvr, work, &lwork,        &info);
	
	// Intialise work space
	lwork = (int) (cplx(work[0]).real());
	free (work); work = (T*) malloc (lwork*sizeof(T));

	// Actual eigen value comp
	if (typeid(T) == typeid(cplx))
		cgeev_ (&jobvl, &jobvr, &n, &this->At(0), &lda, &ev->At(0),     &lev->At(0), &ldvl, &rev->At(0), &ldvr, work, &lwork, rwork, &info);
	else if (typeid(T) == typeid(double))
		dgeev_ (&jobvl, &jobvr, &n, &this->At(0), &lda, &ev->At(0), wi, &lev->At(0), &ldvl, &rev->At(0), &ldvr, work, &lwork,        &info);
	
	// Clean up
	free (wi);
	free (work);
	free (rwork);

	return info;

}


template<class T>
int Matrix<T>::SVD (Matrix<T>* u, Matrix<T>* v, Matrix<T>* s, const char jobz) {

	// SVD only defined on 2D data
	if (!Is2D())
		return -2;
	
	bool    cm   = true;

	int    m     = _dim[COL];
	int    n     = _dim[LIN];
	int    lwork = -1;
	int    info  = 0;
	int    lda   = _dim[COL];
	int    ldu   = (cm) ? ((m >= n) ? m : n) : 1;
	int    ldvt  = (cm) ? ((m >= n) ? n : m) : 1;
	
	T*     work  =     (T*) malloc (sizeof(T));
	float* rwork = (float*) malloc (sizeof(float));
	int*   iwork =   (int*) malloc (8 * MIN(_dim[0],_dim[1]) * sizeof(int));
	
	// Only needed for complex data
	if (typeid(T) == typeid(cplx)) {
		free (rwork);
		rwork = (float*) malloc (MIN(m,n) * MAX(5*MIN(m,n)+7,2*MAX(m,n)+2*MIN(m,n)+1) * sizeof(float));
	}
	
	// Workspace query
	if (typeid(T) == typeid(cplx))
		cgesdd_ (&jobz, &m, &n, &this->At(0), &lda, (float*)&s->At(0), &u->At(0), &ldu, &v->At(0), &ldvt, work, &lwork, rwork, iwork, &info);
	else if (typeid(T) == typeid(double))
		dgesdd_ (&jobz, &m, &n, &this->At(0), &lda,         &s->At(0), &u->At(0), &ldu, &v->At(0), &ldvt, work, &lwork,        iwork, &info);
	
	// Resize work according to ws query
	lwork = (int) cplx (work[0]).real();
	free(work); work = (T*) malloc (lwork * sizeof(T));
	
	//SVD
	if (typeid(T) == typeid(cplx))
		cgesdd_ (&jobz, &m, &n, &this->At(0), &lda, (float*)&s->At(0), &u->At(0), &ldu, &v->At(0), &ldvt, work, &lwork, rwork, iwork, &info);
	else if (typeid(T) == typeid(double))
		dgesdd_ (&jobz, &m, &n, &this->At(0), &lda,         &s->At(0), &u->At(0), &ldu, &v->At(0), &ldvt, work, &lwork,        iwork, &info);
	
	// Clean up
	free (work);
	free (rwork);
	free (iwork);
	
	return info;

} 


template <class T>
Matrix<T>  Matrix<T>::Inv () const {

	Matrix<T> tmp = (*this);

	assert (_dim[COL] == _dim[LIN]);

	int  n    = _dim[COL];

	int info  = 0;

	int* ipiv = (int*) malloc (n * sizeof(int));
	
	char up = 'U';
	char lo = 'L';

	// LQ Factorisation -------------------

	if (typeid (T) == typeid (cplx)) 
		cgetrf_ (&n, &n, &tmp[0], &n, ipiv, &info);
	else if (typeid (T) == typeid (double))
		dgetrf_ (&n, &n, &tmp[0], &n, ipiv, &info);

	// ------------------------------------

	int lwork = -1; 
	T*  work  = (T*) malloc (sizeof(T));

	// Workspace determination ------------

	if (typeid (T) == typeid (cplx)) 
		cgetri_ (&n, &tmp[0], &n, ipiv, work, &lwork, &info);
	else if (typeid (T) == typeid (double))
		dgetri_ (&n, &tmp[0], &n, ipiv, work, &lwork, &info);

	// ------------------------------------

	lwork = (int) cplx (work[0]).real();
	work  = (T*)  realloc (work, lwork * sizeof(T));
	
	// Copy triangular matrix over --------

	if (typeid (T) == typeid (cplx)) 
		cgetri_ (&n, &tmp[0], &n, ipiv, work, &lwork, &info);
	else if (typeid (T) == typeid (double))
		dgetri_ (&n, &tmp[0], &n, ipiv, work, &lwork, &info);

	// ------------------------------------
	
	free (ipiv);
	free (work);

	return tmp;

} 


template<class T>
Matrix<T> 
Matrix<T>::Pinv () {
	
	char      trans  = 'N';

    int       m      = Dim(0);
    int       n      = Dim(1);
    int       nrhs   = m;
	int       lda    = m;
	int       ldb    = MAX(m,n);
	int       lwork  = -1; 
	int       rank   =  0;
	int       info   =  0;
	
	T*        work   = (T*)     malloc (         sizeof(T)); 
	int*      iwork  = (int*)   malloc (         sizeof(int)); 
	float*    rwork  = (float*) malloc (         sizeof(float));    
	float*    s      = (float*) malloc (MIN(m,n)*sizeof(float));    

	Matrix<T> b      =  Matrix<T>::Id(ldb);

	float     rcond  = -1.0;

	if (typeid(T) == typeid(cplx))
		cgelsd_ (&m, &n, &nrhs, &At(0), &lda, &b.At(0), &ldb, s, &rcond, &rank, work, &lwork, rwork, iwork, &info);
	else if (typeid(T) == typeid(double))
		dgelsd_ (&m, &n, &nrhs, &At(0), &lda, &b.At(0), &ldb, s, &rcond, &rank, work, &lwork,        iwork, &info);

	if (typeid(T) == typeid(cplx))
		rwork = (float*) realloc (rwork, (int)rwork[0]*sizeof(float));

	iwork = (int*) realloc (iwork, iwork[0]*sizeof(int));
	lwork = (int) cplx (work[0]).real();
	work  = (T*) realloc (work, lwork * sizeof(T));

	if (typeid(T) == typeid(cplx))
		cgelsd_ (&m, &n, &nrhs, &At(0), &lda, &b.At(0), &ldb, s, &rcond, &rank, work, &lwork, rwork, iwork, &info);
	else if (typeid(T) == typeid(double))
		dgelsd_ (&m, &n, &nrhs, &At(0), &lda, &b.At(0), &ldb, s, &rcond, &rank, work, &lwork,        iwork, &info);

	free (s);
	free (work);
	free (rwork);
	free (iwork);

	return b;
	
}



template<class T>
Matrix<T> 
Matrix<T>::Cholesky (const char uplo) {

	Matrix<T> tmp = (*this);
	int       info = 0;

	dpotrf_ ( &uplo, &tmp.Width(), &tmp[0], &tmp.Height(), &info);
	
	return tmp;

}
