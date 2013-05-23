#ifndef __PIO_HPP__
#define __PIO_HPP__

#include "Matrix.hpp"

inline static void
print (const Matrix<float,MPI>& M, const std::string& name = "pmat", int nout = 6) {
    
	int    m     = M.GHeight();
	int    n     = M.GWidth();
	int    len   = name.length();
    
	float* iwork = (float*) malloc (M.Height() * sizeof(float));
    
	ScalapackTraits<float>::pxlaprnt (&m, &n, M.Data(), &ione, &ione, M.Desc(), &izero, &izero, name.c_str(), &nout, iwork, len);

	free (iwork);
    
}

#endif //__PIO_HPP__
