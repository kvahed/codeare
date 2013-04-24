template<>
struct SSETraits< std::complex<double> > {

    typedef std::complex<double> type;
    typedef __m256d Register;         /**< @brief register type */
    static const unsigned int ne = 2; /**< @brief # of processed elements */
    static const unsigned int ns = 2; /**< @brief # of sub elements */

    /**
     * @brief     AVX load packed aligned
     */
    static inline Register 
    loada (const type* p) {
        return _mm256_load_pd ((double*)p); 
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadu (const type* p) {
        return _mm256_loadu_pd ((double*)p); 
    }

    /**
     * @brief     AVX load packed aligned
     */
    static inline Register 
    loadoa (const type* p) {
        return _mm256_load_pd ((double*)p); 
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadou (const type* p) {
        return _mm256_loadu_pd ((double*)p); 
    }

    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (type* p, Register a) {
        _mm256_store_pd ((double*)p, a); 
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (type* p, Register a) {
        _mm256_storeu_pd ((double*)p, a); 
    }

    /**
     * @brief     AVX packed addition
     */
    static inline Register
    addp (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }

    /**
     * @brief     AVX single addition
     */
    static inline Register 
    adds (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }

    /**
     * @brief     AVX packed subtraction
     */
    static inline Register 
    subp (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }

    /**
     * @brief     AVX single subtraction
     */
    static inline Register 
    subs (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }

    /**
     * @brief     SSE3 packed multiplication
     */
    static inline Register 
    mulp (const Register &a, const Register &b) {

        Register c, d, e, f;

        c = _mm256_mul_pd (a,b);                
        d = _mm256_shuffle_pd (c,c,0x1);        
        e = _mm256_shuffle_pd (b,b,0x1);        
        f = _mm256_sub_pd (c,d);                
        d = _mm256_mul_pd (a,e);                
        e = _mm256_shuffle_pd (d,d,0x1);        
        e = _mm256_add_pd (d,e);    

        return _mm256_shuffle_pd (f,e,0x2);

    }

    /**
     * @brief     AVX single multiplication
     */
    static inline Register 
    muls (const Register &a, const Register &b) {
        return mulp(a, b);
    }

    /**
     * @brief     AVX packed division
     */
    static inline Register 
    divp (const Register &a, const Register &b) {
        return _mm256_div_pd(a, b);
    }

    /**
     * @brief     AVX single division
     */
    static inline Register 
    divs (const Register &a, const Register &b) {
        return _mm256_div_pd(a, b);
    }

    /**
     * @brief     AVX packed SQRT
     */
    static inline Register 
    sqrtp (const Register &a) {
        return _mm256_sqrt_pd(a);
    }

    /**
     * @brief     AVX single SQRT
     */
    static inline Register 
    sqrts (const Register &a, const Register &b) {
        return _mm256_sqrt_pd(a);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    minp (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    mins (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    maxp (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    maxs (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }

}; // SSETraits< std::complex<double> >

template<>
struct SSETraits<double> {

    typedef __m256d Register;         /**< @brief register type */
    static const unsigned int ne = 4; /**< @brief # of processed elements */
    static const unsigned int ns = 1; /**< @brief # of sub elements */

    /**
     * @brief     AVX load packed aligned
     */
    static inline Register 
    loada (const double* p) {
        return _mm256_load_pd (p); 
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadu (const double* p) {
        return _mm256_loadu_pd (p); 
    }

    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (double* p, Register a) {
        _mm256_store_pd (p, a); 
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (double* p, Register a) {
        _mm256_storeu_pd (p, a); 
    }

    /**
     * @brief     AVX packed addition
     */
    static inline Register
    addp (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }

    /**
     * @brief     AVX single addition
     */
    static inline Register 
    adds (const Register &a, const Register &b) {
        return _mm256_add_pd(a, b);
    }

    /**
     * @brief     AVX packed subtraction
     */
    static inline Register 
    subp (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }

    /**
     * @brief     AVX single subtraction
     */
    static inline Register 
    subs (const Register &a, const Register &b) {
        return _mm256_sub_pd(a, b);
    }

    /**
     * @brief     AVX packed multiplication
     */
    static inline Register 
    mulp (const Register &a, const Register &b) {
        return _mm256_mul_pd(a,b);                
    }

    /**
     * @brief     AVX single multiplication
     */
    static inline Register 
    muls (const Register &a, const Register &b) {
        return _mm256_mul_pd(a,b);                
    }

    /**
     * @brief     AVX packed division
     */
    static inline Register 
    divp (const Register &a, const Register &b) {
        return _mm256_div_pd(a, b);
    }

    /**
     * @brief     AVX single division
     */
    static inline Register 
    divs (const Register &a, const Register &b) {
        return _mm256_div_pd(a, b);
    }

    /**
     * @brief     AVX packed SQRT
     */
    static inline Register 
    sqrtp (const Register &a) {
        return _mm256_sqrt_pd(a);
    }

    /**
     * @brief     AVX single SQRT
     */
    static inline Register 
    sqrts (const Register &a, const Register &b) {
        return _mm256_sqrt_pd(a);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    minp (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    mins (const Register &a, const Register &b) {
        return _mm256_min_pd(a, b);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    maxp (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    maxs (const Register &a, const Register &b) {
        return _mm256_max_pd(a, b);
    }

}; // SSETraits< double >

template<>
struct SSETraits<float> {

    typedef float type;
    typedef __m256 Register;         /**< @brief register type */
    static const unsigned int ne = 8; /**< @brief # of processed elements */
    static const unsigned int ns = 1; /**< @brief # of processed elements */

    /**
     * @brief     AVX load packed aligned
     */
    static inline Register 
    loada (const type* p) {
        return _mm256_load_ps (p); 
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadu (const type* p) {
        return _mm256_loadu_ps (p); 
    }

    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (type* p, Register a) {
        _mm256_store_ps (p, a); 
    }

    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (type* p, Register a) {
        _mm256_storeu_ps (p, a); 
    }

    /**
     * @brief     AVX packed addition
     */
    static inline Register
    addp (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }

    /**
     * @brief     AVX single addition
     */
    static inline Register 
    adds (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }

    /**
     * @brief     AVX packed subtraction
     */
    static inline Register 
    subp (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }

    /**
     * @brief     AVX single subtraction
     */
    static inline Register 
    subs (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }

    /**
     * @brief     AVX packed multiplication
     */
    static inline Register 
    mulp (const Register &a, const Register &b) {
        return _mm256_mul_ps(a, b);
    }

    /**
     * @brief     AVX single multiplication
     */
    static inline Register 
    muls (const Register &a, const Register &b) {
        return _mm256_mul_ps(a, b);
    }

    /**
     * @brief     AVX packed division
     */
    static inline Register 
    divp (const Register &a, const Register &b) {
        return _mm256_div_ps(a, b);
    }

    /**
     * @brief     AVX single division
     */
    static inline Register 
    divs (const Register &a, const Register &b) {
        return _mm256_div_ps(a, b);
    }

    /**
     * @brief     AVX packed SQRT
     */
    static inline Register 
    sqrtp (const Register &a) {
        return _mm256_sqrt_ps(a);
    }

    /**
     * @brief     AVX single SQRT
     */
    static inline Register 
    sqrts (const Register &a) {
        return _mm256_sqrt_ps(a);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    minp (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    mins (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    maxp (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    maxs (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }

}; // SSETraits<float>    

template<>
struct SSETraits< std::complex<float> > {
    
    typedef std::complex<float> type;
    typedef __m256 Register;         /**< @brief register type */
    static const unsigned int ne = 4; /**< @brief # of processed elements */
    static const unsigned int ns = 2; /**< @brief # of processed elements */
    
    /**
     * @brief     AVX load packed aligned
     */
    static inline Register 
    loada (const type* p) {
        return _mm256_load_ps ((float*)p); 
    }
    
    /**
     * @brief     AVX load packed unaligned
     */
    static inline Register
    loadu (const type* p) {
        return _mm256_loadu_ps ((float*)p); 
    }
    
    /**
     * @brief     AVX load packed aligned
     */
    static inline void
    stora (type* p, Register a) {
        _mm256_store_ps ((float*)p, a); 
    }
    
    /**
     * @brief     AVX load packed unaligned
     */
    static inline void
    storu (type* p, Register a) {
        _mm256_storeu_ps ((float*)p, a); 
    }
    
    /**
     * @brief     AVX packed addition
     */
    static inline Register
    addp (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }

    /**
     * @brief     AVX single addition
     */
    static inline Register 
    adds (const Register &a, const Register &b) {
        return _mm256_add_ps(a, b);
    }

    /**
     * @brief     AVX packed subtraction
     */
    static inline Register 
    subp (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }

    /**
     * @brief     AVX single subtraction
     */
    static inline Register 
    subs (const Register &a, const Register &b) {
        return _mm256_sub_ps(a, b);
    }

    /**
     * @brief     AVX packed multiplication
     */
    static inline Register 
    mulp (const Register &ab, const Register &cd) {

        Register aa, bb, dc, x0, x1;  

        aa = _mm256_moveldup_ps(ab);  
        bb = _mm256_movehdup_ps(ab);  
        x0 = _mm256_mul_ps(aa, cd);    //ac ad  
        dc = _mm256_shuffle_ps(cd, cd, _MM_SHUFFLE(2,3,0,1));  
        x1 = _mm256_mul_ps(bb, dc);    //bd bc  

        return _mm256_addsub_ps(x0, x1);

    }

    /**
     * @brief     AVX single multiplication
     */
    static inline Register 
    muls (const Register &a, const Register &b) {
        return mulp(a, b);
    }

    /**
     * @brief     AVX packed division
     */
    static inline Register 
    divp (const Register &ab, const Register &cd) {

        Register aa, bb, dc, x0, x1, x2, x3;  

        bb = _mm256_movehdup_ps (ab);
        aa = _mm256_moveldup_ps (ab);
        dc = _mm256_shuffle_ps(cd, cd, _MM_SHUFFLE(2,3,0,1));  
        x0 = _mm256_mul_ps (bb, cd);
        x1 = _mm256_mul_ps (aa, dc);
        x2 = _mm256_addsub_ps (x0, x1);
        x1 = _mm256_mul_ps (dc, dc);
        x0 = _mm256_shuffle_ps (x1, x1, _MM_SHUFFLE(2,3,0,1));
        x3 = _mm256_add_ps (x1, x0);
        x1 = _mm256_div_ps (x2, x3);

        return _mm256_shuffle_ps (x1, x1, _MM_SHUFFLE(2,3,0,1));

    }

    /**
     * @brief     AVX single division
     */
    static inline Register 
    divs (const Register &a, const Register &b) {
        return divp(a, b);
    }

    /**
     * @brief     AVX packed SQRT
     */
    static inline Register 
    sqrtp (const Register &a) {
        return _mm256_sqrt_ps(a);
    }

    /**
     * @brief     AVX single SQRT
     */
    static inline Register 
    sqrts (const Register &a) {
        return _mm256_sqrt_ps(a);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    minp (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    mins (const Register &a, const Register &b) {
        return _mm256_min_ps(a, b);
    }

    /**
     * @brief     AVX packed comparison
     */
    static inline Register 
    maxp (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }

    /**
     * @brief     AVX single comparison
     */
    static inline Register 
    maxs (const Register &a, const Register &b) {
        return _mm256_max_ps(a, b);
    }

}; // SSETraits<float>    
