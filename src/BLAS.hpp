extern "C" {

	// This is all defined only on 2D
	
	// Matrix vector multiplication
	void   dgemv_  (char *trans, int    *m, int *n, void *alpha, void *a, int *lda, void *x, int *incx, void *beta, 
				    void     *y, int *incy);
	void   cgemv_  (char *trans, int    *m, int *n, void *alpha, void *a, int *lda, void *x, int *incx, void *beta, 
				    void     *y, int *incy);
	
	// Matrix matrix multiplication
	void dgemm_    (char *transa, char *transb, int  *m, int   *n, int *k, void *alpha, void *a, int *lda, void *b, 
				    int     *ldb, void   *beta, void *c, int *ldc);
	void cgemm_    (char *transa, char *transb, int  *m, int   *n, int *k, void *alpha, void *a, int *lda, void *b, 
				    int     *ldb, void   *beta, void *c, int *ldc);
	
}

