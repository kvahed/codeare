#ifndef __PIO_HPP__
#define __PIO_HPP__

#include "PMatrix.hpp"

template <class T> inline static void
print (const PMatrix<T>& PM, const std::string& name = "pmat", int nout = 6) {
    
	int   m     = PM.g_m();
	int   n     = PM.g_n();
	int   len   = name.length();

	T*    iwork = (T*) malloc (PM.Height() * sizeof(T));

	ScalapackTraits<T>::pxlaprnt (&m, &n, PM.Data(), &ione, &ione, PM.Desc(), &izero, &izero, name.c_str(), &nout, iwork, len);

	free (iwork);
    
}

#endif //__PIO_HPP__
