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

extern "C" {

	// This is all defined only on 2D
	
	// Matrix vector multiplication
	void   dgemv_  (char *trans, int    *m, int *n, void *alpha, void *a, int *lda, void *x, int *incx, void *beta, void *y, int *incy);
	void   cgemv_  (char *trans, int    *m, int *n, void *alpha, void *a, int *lda, void *x, int *incx, void *beta, void *y, int *incy);
	
	// Matrix matrix multiplication
	void dgemm_    (char *transa, char *transb, int  *m, int   *n, int *k, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
	void cgemm_    (char *transa, char *transb, int  *m, int   *n, int *k, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
	
}

