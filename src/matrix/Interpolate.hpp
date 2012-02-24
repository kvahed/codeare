#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
     
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
interp1 (const Matrix<T>& x, const Matrix<T>& y, const Matrix<double>& xi, const INTERP::Method& intm) {

	Matrix<T> yi;

	gsl_interp_accel *acc 
		= gsl_interp_accel_alloc ();
	gsl_spline *spline 
		= gsl_spline_alloc (gsl_interp_cspline, 10);
	
	gsl_spline_init (spline, x, y, 10);
    
	for (xi = x[0]; xi < x[9]; xi += 0.01) {
		yi = gsl_spline_eval (spline, xi, acc);
		printf ("%g %g\n", xi, yi);
	}

	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
	
	return yi;
	
} 
