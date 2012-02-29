#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "Algos.hpp"
     
namespace INTERP {

	enum Method {

		LINEAR,
		POLYNOMIAL,
		CSPLINE,
		CSPLINE_PERIODIC,
		AKIMA,
		AKIMA_PERIODIC
		
	};

}




template <class T> Matrix<T>
interp1 (Matrix<double>& x, Matrix<T>& y, const Matrix<double>& xi, const INTERP::Method& intm) {

	size_t nx = size(x,0);
	size_t nxi = size(xi,0);

	Matrix<T> yi (nxi,1);

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline;

	if      (intm == INTERP::LINEAR)
		spline = gsl_spline_alloc (gsl_interp_linear, nx);
	else if (intm == INTERP::POLYNOMIAL)
		spline = gsl_spline_alloc (gsl_interp_polynomial, nx);
	else if (intm == INTERP::CSPLINE)
		spline = gsl_spline_alloc (gsl_interp_cspline, nx);
	else if (intm == INTERP::CSPLINE_PERIODIC)
		spline = gsl_spline_alloc (gsl_interp_cspline_periodic, nx);
	else if (intm == INTERP::AKIMA)
		spline = gsl_spline_alloc (gsl_interp_akima, nx);
	else if (intm == INTERP::AKIMA_PERIODIC)
		spline = gsl_spline_alloc (gsl_interp_akima_periodic, nx);
	
	gsl_spline_init (spline, &x[0], &y[0], nx);
    
	for (size_t i = 0; i < xi.Size(); i++)
		yi[i] = gsl_spline_eval (spline, xi[i], acc);

	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
	
	return yi;
	
} 
