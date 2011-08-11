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
int  Matrix<T>::EIG (const bool cv, Matrix<raw>* ev, Matrix<T>* lev, Matrix<T>* rev) {

	// 2D square matrix.

	if (!Is2D())
		return -2;

	if (_dim[COL] != _dim[LIN])
		return -3;
	
#ifdef HAVE_LAPACK

	char    jobvl = (cv) ? 'V' : 'N';
	char    jobvr = (cv) ? 'V' : 'N';
	
	int     n     = _dim[COL];
	
	T*      a     = new T[Size()];
	
	int     i = 0, j = 0, k = 0;

	for (j = 0; j < _dim[COL]; j++)
		for (k = 0; k < _dim[LIN]; k++, i++)
			a[i] = _M[j+k*_dim[LIN]];

	int     lda   = n;
	
	T*       w    = new T[n];

	T*       wi;
	if (typeid(T) == typeid(double))
		wi        = new T[n];
	else 
		wi        = new T[1];

    int     ldvl  =  (cv) ? n : 1;
    int     ldvr  =  (cv) ? n : 1;
	
	T*       vl   = new T[ldvl*ldvl];
	T*       vr   = new T[ldvr*ldvr];

	T*      work  = new T[1];

	int     info  =  0;
	int     lwork = -1;

	float*  rwork;
	if (typeid(T) == typeid(raw))
		rwork = new float[2*n];
	else 
		rwork = new float[1];

	if (typeid(T) == typeid(raw))
		cgeev_ (&jobvl, &jobvr, &n, a, &lda, w,     vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
	else if (typeid(T) == typeid(double))
		dgeev_ (&jobvl, &jobvr, &n, a, &lda, w, wi, vl, &ldvl, vr, &ldvr, work, &lwork,        &info);
	
	lwork = (int) (raw(work[0]).real());

	//std::cout << "xgeev_ (lwork = " << lwork << ")" << std::endl;

	delete [] work;
	work = new T[lwork];

	if (typeid(T) == typeid(raw))
		cgeev_ (&jobvl, &jobvr, &n, a, &lda, w,     vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
	else if (typeid(T) == typeid(double))
		dgeev_ (&jobvl, &jobvr, &n, a, &lda, w, wi, vl, &ldvl, vr, &ldvr, work, &lwork,        &info);

	ev->Dim(COL) = n;
	ev->Reset();
	
	if (typeid(T) == typeid(raw))
		for (i = 0; i < n; i++)
			ev->At(i) = w[i];

	else if (typeid(T) == typeid(double))
		for (i = 0; i < n; i++)
			ev->At(i) = raw(raw(w[i]).real(),raw(wi[i]).real());

	if (cv) {

		lev->Dim(0) = ldvl;
		lev->Dim(1) = ldvl;
		lev->Reset();
		
		i = 0;
		for (j = 0; j < ldvl; j++)
			for (k = 0; k < ldvl; k++, i++)
				lev->At(i) = vl[j+k*ldvl];
		
		rev->Dim(0) = ldvr;
		rev->Dim(1) = ldvr;
		rev->Reset();

		i = 0;
		for (j = 0; j < ldvr; j++)
			for (k = 0; k < ldvr; k++, i++)
				rev->At(i) = vl[j+k*ldvr];
		
	}
		
	delete [] a;
	delete [] w;
	delete [] wi;
	delete [] vl;
	delete [] vr;
	delete [] work;
	delete [] rwork;

	return info;

#else

	return -1;

#endif

}


template<class T>
int Matrix<T>::SVD (const bool cm, Matrix<T>* u, Matrix<T>* v, Matrix<T>* s) {

	if (!Is2D())
		return -2;

#ifndef HAVE_LAPACK

	char    jobz = (cm) ? 'A' : 'N';
	
	int     m    = _dim[COL];
	int     n    = _dim[LIN];
	int    lwork = -1;
	int     info = 0;


	int     lda  = _dim[COL];
	int     ldu  = (cm) ? ((m >= n) ? m : n) : 1;
	int     ldvt = (cm) ? ((m >= n) ? n : m) : 1;

	int i = 0, j = 0, k = 0;

	T*      work =     (T*) malloc (sizeof(T));
	float* rwork = (float*) malloc (sizeof(float));
	int*   iwork =   (int*) malloc (8 * MIN(_dim[0],_dim[1]) * sizeof(int));
	
	if (typeid(T) == typeid(raw)) {
		free (rwork);
		rwork = (float*) malloc (MIN(m,n) * MAX(5*MIN(m,n)+7,2*MAX(m,n)+2*MIN(m,n)+1) * sizeof(float));
	}
	if (typeid(T) == typeid(raw))
		cgesdd_ (&jobz, &m, &n, &this->At(0), &lda, (float*)&s->At(0), &u->At(0), &ldu, &v->At(0), &ldvt, work, &lwork, rwork, iwork, &info);
	else if (typeid(T) == typeid(double))
		dgesdd_ (&jobz, &m, &n, &this->At(0), &lda,         &s->At(0), &u->At(0), &ldu, &v->At(0), &ldvt, work, &lwork,        iwork, &info);
	
	lwork = (int) raw (work[0]).real();
#ifdef DEBUG
	printf ("xgesdd_ (lwork = %i, info = 0)\n", lwork, info) ;
#endif
	free(work); work = (T*) malloc (lwork * sizeof(T));
	
	if (typeid(T) == typeid(raw))
		cgesdd_ (&jobz, &m, &n, &this->At(0), &lda, (float*)&s->At(0), &u->At(0), &ldu, &v->At(0), &ldvt, work, &lwork, rwork, iwork, &info);
	else if (typeid(T) == typeid(double))
		dgesdd_ (&jobz, &m, &n, &this->At(0), &lda,         &s->At(0), &u->At(0), &ldu, &v->At(0), &ldvt, work, &lwork,        iwork, &info);

	free (work);
	free (rwork);
	free (iwork);
	
	return info;

#else //  HAVE_LAPACK
	
	return -1;

#endif // HAVE_LAPACK
	
} 


template <class T>
int  Matrix<T>::Inv () const {

	assert (_dim[COL] == _dim[LIN]);

	int  n    = _dim[COL];

	int info  = 0;

	int* ipiv = (int*) malloc (n * sizeof(int));
	
	char up = 'U';
	char lo = 'L';

	// LQ Factorisation -------------------

	if (typeid (T) == typeid (raw)) 
		cgetrf_ (&n, &n, _M, &n, ipiv, &info);
	else if (typeid (T) == typeid (double))
		dgetrf_ (&n, &n, _M, &n, ipiv, &info);

	// ------------------------------------

	if (info != 0)
		return info;

	int lwork = -1; 
	T*  work  = (T*) malloc (sizeof(T));

	// Workspace determination ------------

	if (typeid (T) == typeid (raw)) 
		cgetri_ (&n, _M, &n, ipiv, work, &lwork, &info);
	else if (typeid (T) == typeid (double))
		dgetri_ (&n, _M, &n, ipiv, work, &lwork, &info);

	// ------------------------------------

	lwork = (int) raw (work[0]).real();
	work  = (T*)  realloc (work, lwork * sizeof(T));
	
	// Copy triangular matrix over --------

	if (typeid (T) == typeid (raw)) 
		cgetri_ (&n, _M, &n, ipiv, work, &lwork, &info);
	else if (typeid (T) == typeid (double))
		dgetri_ (&n, _M, &n, ipiv, work, &lwork, &info);

	// ------------------------------------
	
	free (ipiv);
	free (work);

	return info;

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
	
	T*        work   = (T*) malloc (sizeof(float)); 

	float*    rwork  = (float*) malloc (5*MIN(m,n)*sizeof(float));    
	float*    s      = (float*) malloc (  MIN(m,n)*sizeof(float));    

	Matrix<T> b      =  Matrix<T>::Id(ldb);

	float     rcond  = -1.0;

	if (typeid(T) == typeid(raw))
		cgelss_ (&m, &n, &nrhs, &At(0), &lda, &b.At(0), &ldb, s, &rcond, &rank, work, &lwork, rwork, &info);
	else if (typeid(T) == typeid(double))
		dgelss_ (&m, &n, &nrhs, &At(0), &lda, &b.At(0), &ldb, s, &rcond, &rank, work, &lwork, &info);

	lwork = (int) raw (work[0]).real();
	work  = (T*)  realloc (work, lwork * sizeof(T));

	if (typeid(T) == typeid(raw))
		cgelss_(&m, &n, &nrhs, &At(0), &lda, &b.At(0), &ldb, s, &rcond, &rank, work, &lwork, rwork, &info);
	else if (typeid(T) == typeid(double))
		dgelss_(&m, &n, &nrhs, &At(0), &lda, &b.At(0), &ldb, s, &rcond, &rank, work, &lwork, &info);

	free (rwork);
	free (s);

	return b;
	
}
