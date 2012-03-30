#ifndef __CREATORS_HPP__
#define __CREATORS_HPP__

#include "Matrix.hpp"

#include <gsl/gsl_rng.h>
#include <limits>


/**
 * @brief       Zero matrix
 *
 * @param  col  Column
 * @param  lin  Rows
 * @param  cha  Dimension
 * @param  set  Dimension
 * @param  eco  Dimension
 * @param  phs  Dimension
 * @param  rep  Dimension
 * @param  seg  Dimension
 * @param  par  Dimension
 * @param  slc  Dimension
 * @param  ida  Dimension
 * @param  idb  Dimension
 * @param  idc  Dimension
 * @param  idd  Dimension
 * @param  ide  Dimension
 * @param  ave  Dimension
 *
 * @return      Zero matrix
 *
 */
template <class T> inline static Matrix<T> 
zeros           (const size_t& col, 
				 const size_t& lin = 1, 
				 const size_t& cha = 1,
				 const size_t& set = 1,
				 const size_t& eco = 1,
				 const size_t& phs = 1,
				 const size_t& rep = 1,
				 const size_t& seg = 1,
				 const size_t& par = 1,
				 const size_t& slc = 1,
				 const size_t& ida = 1,
				 const size_t& idb = 1,
				 const size_t& idc = 1,
				 const size_t& idd = 1,
				 const size_t& ide = 1,
				 const size_t& ave = 1) {

 	return Matrix<T> (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave);

}


/**
 * @brief       Zero matrix
 *
 * @param  sz   Size vector
 * @return      Zero matrix
 *
 */
template <class T> inline static Matrix<T> 
zeros           (const Matrix<size_t>& sz) {

	size_t n [INVALID_DIM], i;

	for (i = 0; i < numel(sz) && i < INVALID_DIM; i++)
		n[i] = sz[i];

	for (     ; i < INVALID_DIM; i++)
		n[i] = 1; 

 	return Matrix<T> (n);

}



/**
 * @brief      Ones matrix
 *
 * @param  col  Column
 * @param  lin  Rows
 * @param  cha  Dimension
 * @param  set  Dimension
 * @param  eco  Dimension
 * @param  phs  Dimension
 * @param  rep  Dimension
 * @param  seg  Dimension
 * @param  par  Dimension
 * @param  slc  Dimension
 * @param  ida  Dimension
 * @param  idb  Dimension
 * @param  idc  Dimension
 * @param  idd  Dimension
 * @param  ide  Dimension
 * @param  ave  Dimension
 *
 * @return      Ones matrix
 *
 */
template <class T> inline static Matrix<T> 
ones            (const size_t& col, 
				 const size_t& lin = 1, 
				 const size_t& cha = 1,
				 const size_t& set = 1,
				 const size_t& eco = 1,
				 const size_t& phs = 1,
				 const size_t& rep = 1,
				 const size_t& seg = 1,
				 const size_t& par = 1,
				 const size_t& slc = 1,
				 const size_t& ida = 1,
				 const size_t& idb = 1,
				 const size_t& idc = 1,
				 const size_t& idd = 1,
				 const size_t& ide = 1,
				 const size_t& ave = 1) {

 	 Matrix<T> res (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave);
	 size_t i = numel(res);

	 while (i--)
		 res[i] = T(1);

	 return res;

}



/**
 * @brief       Random matrix
 *
 * @param  col  Column
 * @param  lin  Rows
 * @param  cha  Dimension
 * @param  set  Dimension
 * @param  eco  Dimension
 * @param  phs  Dimension
 * @param  rep  Dimension
 * @param  seg  Dimension
 * @param  par  Dimension
 * @param  slc  Dimension
 * @param  ida  Dimension
 * @param  idb  Dimension
 * @param  idc  Dimension
 * @param  idd  Dimension
 * @param  ide  Dimension
 * @param  ave  Dimension
 *
 * @return      Random matrix
 *
 */
template<class T> static Matrix<T> 
rand           (const size_t& col, 
				const size_t& lin = 1, 
				const size_t& cha = 1,
				const size_t& set = 1,
				const size_t& eco = 1,
				const size_t& phs = 1,
				const size_t& rep = 1,
				const size_t& seg = 1,
				const size_t& par = 1,
				const size_t& slc = 1,
				const size_t& ida = 1,
				const size_t& idb = 1,
				const size_t& idc = 1,
				const size_t& idd = 1,
				const size_t& ide = 1,
				const size_t& ave = 1) {
	
	Matrix<T> res (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave);

	size_t i = numel(res);
	
	const gsl_rng_type* grt;
	gsl_rng* r;
	
	gsl_rng_env_setup();
	grt = gsl_rng_default;
	r = gsl_rng_alloc (grt);

	if      (typeid(T) == typeid(float) || typeid(T) == typeid(double))
		while (i--)
			res[i] = gsl_rng_uniform (r);
	else if (typeid(T) == typeid(cxdb))
		while (i--) {
			((double*) &res[i])[0] = gsl_rng_uniform (r);
			((double*) &res[i])[1] = gsl_rng_uniform (r);
		}
	else if (typeid(T) == typeid(cxfl))
		while (i--) {
			((float*) &res[i])[0] = gsl_rng_uniform (r);
			((float*) &res[i])[1] = gsl_rng_uniform (r);
		}
	else if (typeid(T) == typeid(short))
		while (i--)
			res[i] = gsl_rng_uniform_int (r, SHRT_MAX);
	
	else if (typeid(T) == typeid(long))
		while (i--)
			res[i] = gsl_rng_uniform_int (r, INT_MAX);
	
	gsl_rng_free (r);

	return res;

}



/**
 * @brief       nxn square matrix with circle centered at p
 *
 * @param  p    Center point of circle
 * @param  n    Side length of square
 * @param  s    Scaling factor
 * @return      Matrix with circle
 */
template <class T> inline static Matrix<T>
circle (const float* p, const size_t n, const T s = T(1)) {

	Matrix<T> res(n);

	float m[2];
	float rad;

	rad = p[0] * float(n) / 2.0;

	m[0] = (1.0 - p[1]) * float(n) / 2.0;
	m[1] = (1.0 - p[2]) * float(n) / 2.0;

	for (size_t r = 0; r < res.Dim(1); r++)
		for (size_t c = 0; c < res.Dim(0); c++)
			res(c,r) = ( pow(((float)c-m[0])/rad, 2.0 ) + pow(((float)r-m[0])/rad, 2.0) <= 1.0) ? s : T(0.0);

	return res;

}



/**
 * @brief       nxnxn cube with sphere centered at p
 *
 * @param  p    Center point of sphere
 * @param  n    Side length of cube
 * @param  s    Scaling factor
 * @return      Matrix with circle
 */
template <class T> inline static Matrix<T>
sphere (const float* p, const size_t n, const T s = T(1)) {

	Matrix<T> res (n,n,n);

	float m[3];
	float rad;

	rad = p[0] * float(n) / 2.0;

	m[0] = (1.0 - p[1]) * float(n) / 2.0;
	m[1] = (1.0 - p[2]) * float(n) / 2.0;
	m[2] = (1.0 - p[3]) * float(n) / 2.0;

	for (size_t s = 0; s < res.Dim(2); s++)
		for (size_t r = 0; r < res.Dim(1); r++)
			for (size_t c = 0; c < res.Dim(0); c++)
				res(c,r) = ( pow (((float)c-m[0])/rad, 2.0) + pow (((float)r-m[1])/rad, 2.0) + pow (((float)s-m[2])/rad, 2.0) <= 1.0) ? s : T(0.0);

	return res;

}



/**
 * @brief       nxn square matrix with circle centered at p
 *
 * @param  p    Center point of ellipse and excentricities
 * @param  n    Side length of square
 * @param  s    Scaling 
 * @return      Matrix with circle
 */
template <class T> inline static Matrix<T>
ellipse (const float* p, const size_t n, const T s = T(1)) {

	Matrix<T> res (n);

	float m[2];
	float a[2];

	a[0] = p[0] * float(n) / 2.0;
	a[1] = p[1] * float(n) / 2.0;

	m[0] = (1.0 - p[2]) * float(n) / 2.0;
	m[1] = (1.0 - p[3]) * float(n) / 2.0;

	float cosp = cos(p[4]);
	float sinp = sin(p[4]);
	
#pragma omp parallel default (shared) 
	{
		
		size_t tid      = omp_get_thread_num();
		size_t chunk    = n / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk) 
		
	for (size_t r = 0; r < n; r++)
		for (size_t c = 0; c < n; c++)
			res(c,r) = (pow( (((float)c-m[1])*cosp+((float)r-m[0])*sinp)/a[1], 2.0 ) + 
						pow( (((float)r-m[0])*cosp-((float)c-m[1])*sinp)/a[0], 2.0) <= 1.0) ? s : T(0.0);

	}

	return res;

}



/**
 * @brief       nxnxn cube with ellipsoid centered at p
 *
 * @param  p    Center point of ellipsoid and excentricities
 * @param  n    Side length of square
 * @param  s    Scaling 
 * @return      Cube with ellipsoid
 */
template <class T> inline static Matrix<T>
ellipsoid (const float* p, const size_t n, const T s) {

	Matrix<T> res (n,n,n);

	float m[3];
	float a[3];
	float d;

	a[0] = p[0] * float(n) / 2.0;
	a[1] = p[1] * float(n) / 2.0;
	a[2] = p[2] * float(n) / 2.0;

	m[0] = (1.0 - p[3]) * float(n) / 2.0;
	m[1] = (1.0 - p[4]) * float(n) / 2.0;
	m[2] = (1.0 - p[5]) * float(n) / 2.0;

	float cosp = cos(p[6]);
	float sinp = sin(p[6]);
	
#pragma omp parallel default (shared) 
	{
		
		size_t tid      = omp_get_thread_num();
		size_t chunk    = n / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk) 
		
		for (size_t s = 0; s < n; s++)
			for (size_t r = 0; r < n; r++)
				for (size_t c = 0; c < n; c++)
					res(c,r,s) = ( pow( (((float)c-m[1])*cosp+((float)r-m[0])*sinp)/a[1], 2.0) + 
								   pow( (((float)r-m[0])*cosp-((float)c-m[1])*sinp)/a[0], 2.0) +
								   pow( ((float)s-m[2])/a[2], 2.0) <= 1.0) ? s : T(0.0);
		
	}

	return res;

}




/**
 * @brief           nxn Shepp-Logan phantom.
 *
 *                  Shepp et al.<br/> 
 *                  The Fourier reconstruction of a head section.<br/> 
 *                  IEEE TNS. 1974; 21: 21-43
 *
 * @param  n        Side length of matrix
 * @return          Shepp-Logan phantom
 */
template<class T> inline static Matrix<T> 
phantom (const size_t& n) {
	
	const size_t ne = 10; // Number of ellipses
	const size_t np = 5;  // Number of geometrical parameters
	
	float p[ne][np] = {
		{ 0.6900, 0.9200,  0.00,  0.0000,  0.0 },
		{ 0.6624, 0.8740,  0.00, -0.0184,  0.0 },
        { 0.1100, 0.3100, -0.22,  0.0000, -0.3 },
		{ 0.1600, 0.4100,  0.22,  0.0000,  0.3 },
		{ 0.2100, 0.2500,  0.00,  0.3500,  0.0 },
		{ 0.0460, 0.0460,  0.00,  0.1000,  0.0 },
		{ 0.0460, 0.0460,  0.00, -0.1000,  0.0 },
		{ 0.0460, 0.0230,  0.08, -0.6050,  0.0 },
		{ 0.0230, 0.0230,  0.00, -0.6060,  0.0 },
		{ 0.0230, 0.0460, -0.06, -0.6050,  0.0 }
	};

	// Size_Tensities
	T v[ne] = {T(1.0), T(-0.8), T(-0.2), T(-0.2), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1)};

	// Empty matrix
	Matrix<T> res (n);
	Matrix<T> e;

	for (size_t i = 0; i < ne; i++) {
		e    = ellipse<T> (p[i], n, v[i]);
		res += e;
	}

	return res;

}





/**
 * @brief           nxnxn Shepp-Logan phantom
 * 
 *                  Koay et al.<br/>
 *                  Three dimensional analytical magnetic resonance imaging phantom in the Fourier domain.<br/>
 *                  MRM. 2007; 58: 430-436
 *
 * @param  n        Side length of matrix
 * @return          nxn zeros
 */
template <class T> inline static Matrix<T> 
phantom3D (const size_t& n) {

	const size_t ne = 10; // Number of ellipses
	const size_t np =  9; // Number of geometrical parameters

	float p[ne][np] = {
		{ 0.690, 0.920, 0.900,  0.00,  0.000,  0.000,  0.0, 0.0, 0.0 },
        { 0.662, 0.874, 0.880,  0.00,  0.000,  0.000,  0.0, 0.0, 0.0 },
        { 0.110, 0.310, 0.220, -0.22,  0.000, -0.250, -0.3, 0.0, 0.0 },
        { 0.160, 0.410, 0.210,  0.22,  0.000, -0.250,  0.3, 0.0, 0.0 },
        { 0.210, 0.250, 0.500,  0.00,  0.350, -0.250,  0.0, 0.0, 0.0 },
        { 0.046, 0.046, 0.046,  0.00,  0.100, -0.250,  0.0, 0.0, 0.0 },
        { 0.046, 0.023, 0.020,  0.08, -0.650, -0.250,  0.0, 0.0, 0.0 },
        { 0.046, 0.023, 0.020,  0.06, -0.650, -0.250,  0.0, 0.0, 0.0 },
        { 0.056, 0.040, 0.100, -0.06, -0.105,  0.625,  0.0, 0.0, 0.0 },
        { 0.056, 0.056, 0.100,  0.00,  0.100,  0.625,  0.0, 0.0, 0.0 }
	};

	T v[ne] = {2.0, -0.8, -0.2, -0.2, 0.2, 0.2, 0.1, 0.1, 0.2, -0.2};

	Matrix<T> res = zeros<T>(n,n,n);
	Matrix<T> e;
	
	for (size_t i = 0; i < ne; i++) {
		e    = ellipsoid<T> (p[i], n, v[i]);
		res += e;
	}

	return res;

}



template<class T> inline static Matrix<size_t>
meshgrid (const Matrix<size_t>& d) {
	
	size_t side [3];
	
	side[0] = d(0,1) - d(0,0) + 1;
	side[1] = d(1,1) - d(1,0) + 1;
	side[2] = d(2,1) - d(2,0) + 1;
	
    Matrix<size_t> mg (side[1], side[0], side[2], 3);
	
	for (size_t s = 0; s < side[2]; s++)
		for (size_t l = 0; l < side[0]; l++)
			for (size_t c = 0; c < side[1]; c++) {
				mg(c,l,s,0) = l + d(0,0);
				mg(c,l,s,1) = c + d(0,1);
				mg(c,l,s,2) = s + d(0,2);
			}
	
	return mg;
	
}



template <class T> inline static Matrix<T>
eye (const size_t n) {

 	Matrix<T> M (n);

 	for (size_t i = 0; i < n; i++)
 		M[i*n+i] = T(1.0);

 	return M;

}



template <class T> inline static Matrix<T> 
linspace (const T& start, const T& end, const size_t& n) {
	
	assert (n > 1);
	
	Matrix<T> res (n, 1);
	T gap;

	gap      = T(end-start) / T(n-1);
	
	res[0]   = start;
	res[n-1] = end;
	
	for (int i = 1; i < n-1; i++)
		res[i] = res[i-1] + gap;
	
	return res;
	
}




#endif


