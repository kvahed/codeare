#ifndef __CREATORS_HPP__
#define __CREATORS_HPP__

#include "Matrix.hpp"
#include "Algos.hpp"
#if !defined(_MSC_VER) || _MSC_VER>1200
#  include "RandTraits.hpp"
#endif



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
		 const size_t& lin,
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
zeros           (const Vector<size_t>& sz) {
 	return Matrix<T> (sz);
}

/**
 * @brief       Square matrix of zeros
 *
 * @param  n    Side length
 * @return      Zero matrix
 */
template <class T> inline static Matrix<T>
zeros            (const size_t& n) {
	return zeros<T>(n,n);
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
				 const size_t& lin,
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
     std::fill (res.Begin(), res.End(), T(1));

	 return res;

}


/**
 * @brief       Square matrix of ones
 *
 * @param  n    Side length
 * @return      Ones matrix
 */
template <class T> inline static Matrix<T>
ones            (const size_t& n) {
	return ones<T>(n,n);
}


/**
 * @brief       Zero matrix
 *
 * @param  sz   Size vector
 * @return      Zero matrix
 *
 */
template <class T> inline static Matrix<T>
ones           (const Vector<size_t>& sz) {
 	return Matrix<T>(sz) = (T)1;
}


#if !defined(_MSC_VER) || _MSC_VER>1200
/**
 * @brief       Uniformly random matrix
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
				const size_t& lin,
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
    Random<T>::Uniform(res);
	return res;

}


/**
 * @brief       Uniformly random matrix
 *
 * @param  sz   Size vector
 * @return      Rand matrix
 *
 */
template <class T> inline static Matrix<T>
rand           (const Vector<size_t>& sz) {

	Matrix<T> res (sz);
    Random<T>::Uniform(res);
 	return res;

}

/**
 * @brief       Random square matrix
 *
 * @param  n    Side length
 * @return      Random matrix
 */
template<class T> static Matrix<T>
rand (const size_t n) {
	return rand<T>(n,n);
}


/**
 * @brief       Uniformly random matrix
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
randn          (const size_t& col, 
		const size_t& lin,
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
	Random<T>::Normal(res);
	return res;

}


/**
 * @brief       Uniformly random matrix
 *
 * @param  sz   Size vector
 * @return      Rand matrix
 *
 */
template <class T> inline static Matrix<T>
randn          (const Vector<size_t>& sz) {

	Matrix<T> res (sz);
	Random<T>::Normal(res);
 	return res;

}

/**
 * @brief       Random square matrix
 *
 * @param  n    Side length
 * @return      Random matrix
 */
template<class T> static Matrix<T>
randn (const size_t n) {
	return rand<T>(n,n);
}
#endif

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
			res(c,r) = ( pow(((float)c-m[0])/rad, (float)2.0 ) + pow(((float)r-m[0])/rad, (float)2.0) <= 1.0) ? s : T(0.0);

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
				res(c,r) = ( pow (((float)c-m[0])/rad, (float)2.0) + pow (((float)r-m[1])/rad, (float)2.0) + pow (((float)s-m[2])/rad, (float)2.0) <= 1.0) ? s : T(0.0);

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
#pragma omp for schedule (dynamic, n / omp_get_num_threads())
		
	for (int r = 0; r < (int)n; r++)
		for (size_t c = 0; c < n; c++) {
			float x = (((float)c-m[1])*cosp+((float)r-m[0])*sinp)/a[1];
			float y = (((float)r-m[0])*cosp-((float)c-m[1])*sinp)/a[0];

			res(c,r) = (x*x + y*y) <= 1.0 ? s : T(0.0);
		}
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
#pragma omp for schedule (dynamic, n / omp_get_num_threads())
		
		for (int s = 0; s < n; s++)
			for (size_t r = 0; r < n; r++)
				for (size_t c = 0; c < n; c++) {
					float x = (((float)c-m[1])*cosp+((float)r-m[0])*sinp)/a[1];
					float y = (((float)r-m[0])*cosp-((float)c-m[1])*sinp)/a[0];
					float z =  ((float)s-m[2])/a[2];
					res(c,r,s) = (x*x + y*y + z*z) <= 1.0 ? s : T(0.0);
				}
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
		{ .69f,   .92f,   .0f,   .0f,     .0f },
		{ .6624f, .874f,  .0f,  -.0184f,  .0f },
        { .11f,   .31f,  -.22f,  .0f,    -.3f },
		{ .16f,   .41f,   .22f,  .0f,     .3f },
		{ .21f,   .25f,   .0f,   .35f,    .0f },
		{ .046f,  .046f,  .00f,  .1f,     .0f },
		{ .046f,  .046f,  .0f,  -.1f,     .0f },
		{ .046f,  .023f,  .08f, -.605f,   .0f },
		{ .023f,  .023f,  .0f,  -.606f,   .0f },
		{ .023f,  .046f, -.06f, -.605f,   .0f }
	};

	// Size_Tensities
#pragma warning (disable : 4305)
	T v[ne] = {(T)1., (T)-.8, (T)-.2, (T)-.2, (T).1, (T).1, (T).1, (T).1, (T).1, (T).1};
#pragma warning (default : 4305)

	// Empty matrix
	Matrix<T> res (n);
	Matrix<T> e;

	for (size_t i = 0; i < ne; i++) {
		e    = ellipse<T> (p[i], n, v[i]);
		res += e;
	}

	return res;

}

template<class T> inline static Matrix<T>
phantom (const size_t& n, const size_t& m) {
	assert (n==m);
	return phantom<T>(n);
}



/**
 * @brief           nxnxn Shepp-Logan phantom
 * 
 *                  Koay et al.<br/>
 *                  Three dimensional analytical magnetic resonance imaging phantom in the Fourier domain.<br/>
 *                  MRM. 2007; 58: 430-436
 *
 * @param  n        Side length o matrix
 * @return          nxn zeros
 */
template <class T> inline static Matrix<T> 
phantom3D (const size_t& n) {

	const size_t ne = 10; // Number o ellipses
	const size_t np =  9; // Number o geometrical parameters

#pragma warning( disable : 4838)
    float p[ne][np] = {
		{ .69,  .92,  .9,   .0,   .0,   .0,   .0, .0, .0 },
        { .662, .874, .88,  .0,   .0,   .0,   .0, .0, .0 },
        { .11,  .31,  .22, -.22,  .0,  -.25, -.3, .0, .0 },
        { .16,  .41,  .21,  .22,  .0,  -.25,  .3, .0, .0 },
        { .21,  .25,  .5,   .0,   .35, -.25,  .0, .0, .0 },
        { .046, .046, .046, .0,   .1,  -.25,  .0, .0, .0 },
        { .046, .023, .02,  .08, -.65, -.25,  .0, .0, .0 },
        { .046, .023, .02,  .06, -.65, -.25,  .0, .0, .0 },
        { .056, .04,  .1,  -.06, -.105, .625, .0, .0, .0 },
        { .056, .056, .1,   .0,   .1,   .625, .0, .0, .0 }
	};
#pragma warning( default : 4838)


	double v[ne] = {2., -.8, -.2, -.2, .2, .2, .1, .1, .2, -.2};

	Matrix<T> res = zeros<T>(n,n,n);
	Matrix<T> e;
	
	for (size_t i = 0; i < ne; i++) {
		e    = ellipsoid<T> (p[i], n, v[i]);
		res += e;
	}

	return res;

}
template<class T> inline static Matrix<T>
phantom (const size_t& n, const size_t& m, const size_t& l) {
	assert (n==m);
	assert (n==l || l==1);
	if (l==n)
		return phantom3D<T>(n);
	else
		return phantom<T>(n);
}



template <class T>
inline static Matrix<T>
eye (const size_t n) {

 	Matrix<T> M (n);

 	for (size_t i = 0; i < n; i++)
 		M[i*n+i] = T(1.0);

 	return M;

}



template <class T> inline static Matrix<T> 
linspace (const T& start, const T& end, const size_t& n) {
	
	assert (n >= 1);
	
	Matrix<T> res (n, 1);
	T gap;

	gap      = T(end-start) / T(n-1);
	
	res[0]   = start;
	res[n-1] = end;
	
	for (size_t i = 1; i < n-1; i++)
		res[i] = res[i-1] + gap;
	
	return res;
	
}


/**
 * @brief    MATLAB-like meshgrid. x and y vectors must be specified z may be specified optionally.
 *
 * @param x  X-Vector
 * @param y  Y-Vector
 * @param z  Z-Vector (default: unused)
 * @return   Mesh grid O (Ny x Nx x Nz x 3) (if z specified) else O (Ny x Nx x 2)<br/>
 */
template <class T> inline static Matrix<T>
meshgrid (const Vector<T>& x, const Vector<T>& y, const Vector<T>& z = Vector<T>(1)) {

	size_t nx = numel(x);
	size_t ny = numel(y);
	size_t nz = numel(z);

	assert (nx > 1);
	assert (ny > 1);

	// Column vectors
	assert (size(x,0) == nx); 
	assert (size(y,0) == ny);
	assert (size(z,0) == nz);

	Matrix<T> res (ny, nx, (nz > 1) ? nz : 2, (nz > 1) ? 3 : 1);
	
	for (size_t i = 0; i < ny * nz; i++) 
		Row    (res, i          , x);
	for (size_t i = 0; i < nx * nz; i++) 
		Column (res, i + nx * nz, y);
	if (nz > 1)
		for (size_t i = 0; i < nz; i++)
			Slice  (res, i +  2 * nz, z[i]);
	
	return res;	

}


template<class T> inline static Matrix<T>
zpad (const Matrix<T>& a, size_t m, size_t n) {
	assert(is2d(a));
	size_t am = size(a,0), an = size(a,1);
	assert(am<=m);
	assert(an<=n);
	size_t am2 = (m-am)/2, an2 = (n-an)/2;
	Matrix<T> ret(m,n);
	for (size_t i = 0; i < an; ++i)
		std::copy(&a(0,i), &a(0,i)+am, &ret(am2,i+an2));
	return ret;
}

template<class T> inline static Matrix<T>
zpad (const Matrix<T>& a, size_t m, size_t n, size_t o) {
	assert(is3d(a));
	size_t am = size(a,0), an = size(a,1), ao = size(a,2);
	assert(am<=m);
	assert(an<=n);
	assert(ao<=o);
	size_t am2 = (m-am)/2, an2 = (n-an)/2, ao2 = (o-ao)/2;
	Matrix<T> ret(m,n,o);
	for (size_t j = 0; j < ao; ++j)
		for (size_t i = 0; i < an; ++i)
			std::copy(&a(0,i,j), &a(0,i,j)+am, &ret(am2,i+an2,j+ao2));
	return ret;
}

template<class T> inline static Matrix<T>
zpad (const Matrix<T>& a, size_t m, size_t n, size_t o, size_t p) {
	assert(is4d(a));
	size_t am = size(a,0), an = size(a,1), ao = size(a,2), ap = size(a,3);
	assert(am<=m);
	assert(an<=n);
	assert(ao<=o);
	assert(ap<=p);
	size_t am2 = (m-am)/2, an2 = (n-an)/2, ao2 = (o-ao)/2, ap2 = (p-ap)/2;
	Matrix<T> ret(m,n,o,p);
	for (size_t k = 0; k < ap; ++k)
		for (size_t j = 0; j < ao; ++j)
			for (size_t i = 0; i < an; ++i)
				std::copy(&a(0,i,j,k), &a(0,i,j,k)+am, &ret(am2,i+an2,j+ao2,k+ap2));
	return ret;
}

template<class T> inline static Matrix<T> repmat (const Matrix<T>& M, const size_t m,
                                                 const size_t n) {
//    assert (is2d(M));
    assert (m>=1);
    assert (n>=1);
    Vector<size_t> odims = size(M), ndims = odims;
    ndims[0] *= m;
    ndims[1] *= n;
    Matrix<T> ret(ndims);
    for (size_t j = 0; j < n*odims[1]; ++j)
        for (size_t i = 0; i < m*odims[0]; ++i)
            ret(i,j) = M(i%odims[0],j%odims[1]);
            //std::copy(M.Begin()+j*odims[0],M.Begin()+(j+1)*odims[0],ret.Begin()+i*odims[0]+j*ndims[0]);
    return ret;
}

#endif


