#include <functional>

/**
     * @name            Some operators
     *                  Operator definitions. Needs big expansion still.
     */
    
    //@{
    


    template<paradigm Q> Matrix<T,P>&
    operator= (const Matrix<T,Q>&);

    
    /**
     * @brief           Assignment data
     *
     * @param  v        Data vector (size must match numel(M)).
     */
    inline Matrix<T,P>&
    operator=           (const container<T>& v) {

    	assert (_M.size() == v.size());

        if (&_M != &v)
            _M = v;

        return *this;

    }



    /**
     * @brief           Assignment operator. Sets all elements s.
     *
     * @param  s        The assigned scalar.
     */
    inline Matrix<T,P>&
    operator=           (const T s) {

        T t = T(s);

#if defined USE_VALARRAY
        _M = t;        
#else
#ifdef EW_OMP
    #pragma omp parallel for
        for (size_t i = 0; i < Size(); ++i)
            _M[i] = t;
#endif
        std::fill(_M.begin(), _M.end(), t);
#endif
        
        return *this;
    }
    
    
    /**
     * @brief           Unary minus (additive inverse)
     *
     * @return          Negation
     */
    inline Matrix<T,P>
    operator-           () const {

        Matrix<T,P> res (_dim);

#if defined USE_VALARRAY
        res = -_M;        
#else
#ifdef EW_OMP
     #pragma omp parallel for
        for (size_t i = 0; i < Size(); ++i)
			res[i] = -_M[i];
#else
        std::transform (_M.begin(), _M.end(), res.Begin(), std::negate<T>());
#endif
#endif
        return res;

    }


    /**
     * @brief           Unary plus
     *
     * @return          Identity
     */
    inline Matrix<T,P>
    operator+           () const {

        return *this;

    }
    
    
    /**
     * @brief           Transposition / Complex conjugation. i.e. this'.
     *
     * @return          Matrix::tr()
     */
    inline Matrix<T,P>
    operator!           () const {

        Matrix<T,P> res (_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < _dim[1]; ++i)
			for (size_t j = 0; j < _dim[0]; j++)
				res(i,j) = this->At(j,i);

        return res;

    }

    
    
    /**
     * @brief           Return a matrix with result[i] = (m[i] ? this[i] : 0).
     *
     * @param  M        The operand
     * @return          Cross-section or zero
     */
    inline Matrix<T,P>
    operator&           (const Matrix<cbool>& M) const ;
    
    
     /**
     * @brief           Scalar inequality. result[i] = (this[i] != m). i.e. this ~= m
     *
     * @param  s        Comparing scalar.
     * @return          Matrix of false where elements are equal s and true else.
     */
    inline Matrix<cbool>
    operator!=          (const T s) const {

        Matrix<cbool> res(_dim);
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = (_M[i] != s) ? 1 : 0;
        return res;

    }

    
    
    /**
     * @brief           Scalar greater comaprison, result[i] = (this[i] > m). i.e. this > m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator>           (const T s) const {

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::greater(_M[i], s);

        return res;

    }

    
    
    /**
     * @brief           Scalar greater or equal comparison. result[i] = (this[i] >= m). i.e. this >= m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator>=          (const T s) const {

		Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::greater_or_equal(_M[i], s);

        return res;

	}

    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] <= m). i.e. this <= m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator<=          (const T s) const {

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::less_or_equal(_M[i], s);

        return res;

    }

    
    /**
     * @brief           Scalar minor or equal comparison. result[i] = (this[i] < m). i.e. this < m
     *
     * @param  s        Comparing scalar.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator<           (const T s) const {

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::less(_M[i], s);

        return res;

    }


    /**
     * @brief           Elementwise equality, result[i] = (this[i] != m[i]). i.e. this ~= m
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator!=          (const Matrix<T,P>& M) const {

        assert (_dim == M._dim);

        Matrix<cbool> res(_dim,_res);
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = (_M[i]!= M[i]) ? 1 : 0;
        return res;

    }

    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] >= m[i]). i.e. this >= m
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator>=          (const Matrix<T,P>& M) const {

        assert (_dim == M._dim);

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::greater_or_equal(_M[i], M[i]);

        return res;

    }

    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] <= m[i]). i.e. this <= m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator<=          (const Matrix<T,P>& M) const {

        assert (_dim == M._dim);

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::less_or_equal(_M[i], M[i]);

        return res;

    }

    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] > m[i]). i.e. this > m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator>           (const Matrix<T,P>& M) const {

        assert (_dim == M._dim);

        Matrix<cbool> res(_dim,_res);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::greater(_M[i], M[i]);

        return res;

    }

    
    /**
     * @brief           Matrix comparison, result[i] = (this[i] < m[i]). i.e. this < m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator<           (const Matrix<T,P>& M) const {

        assert (_dim == M._dim);

        Matrix<cbool> res(_dim,_res);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::less(_M[i], M[i]);

        return res;

    }

    
    /**
     * @brief           Matrix comparison, result[i] = (m[i] || this[i] ? 1 : 0). i.e. this | m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator||          (const Matrix<T,P>& M) const {

        assert (_dim == M._dim);

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::logical_or(_M[i], M[i]);

        return res;

    }

    
    
    /**
     * @brief           Matrix comparison, result[i] = (m[i] && this[i] ? 1 : 0). i.e. this & m.
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator&&          (const Matrix<T,P>& M) const {

        assert (_dim == M._dim);

        Matrix<cbool> res(_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < Size(); ++i)
        	res[i] = CompTraits<T>::logical_and(_M[i], M[i]);

        return res;

    }


    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Power.
     * @return          Result
     */
    inline Matrix<T,P>
    operator^           (const float p) const {

    	Matrix<T,P> res = *this;

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] = (p == 0) ? T(1) : pow(res[i],  p);

        return res;

    }


    /**
     * @brief           Elementwise raise of power. i.e. this .^ p.
     *
     * @param  p        Power.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator^=          (const float p) {

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			_M[i] = pow(_M[i],  p);

        return *this;

    }
    

    /**
     * @brief           Elementwise addition. i.e. this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>
    operator+          (const Matrix<T,P> &M) const {
        
		Matrix<T,P> res = *this;
		return res += M;
        
    }
    
    
    /**
     * @brief           Elementwise addition of two matrices
     *
     * @param  M        Matrix additive.
     */
    template <class S>
    inline const Matrix<T,P>
    operator+          (const Matrix<S,P>& M) const {
        
		Matrix<T,P> res = *this;
		return res += M;
		
    }
    
    
    /**
     * @brief           Elementwise addition iof all elements with a scalar
     *
     * @param  s        Scalar additive.
     */
    template <class S>
    inline Matrix<T,P>
    operator+           (const S s) const {

        Matrix<T,P> res = *this;
        return res += s;

    }

    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator+=         (const Matrix<T,P>& M) {

        assert (_dim == M._dim);

#if defined USE_VALARRAY
        _M += M.Container();        
#elif defined EXPLICIT_SIMD
        if (fp_type(_M[0]))
        	SSE::binary<T>(_M, M.Container(), SSE::add<T>(), _M);
        else
        	for (size_t i = 0; i < Size(); ++i)
        		_M[i] += M[i];
#else
#ifdef EW_OMP
    #pragma omp parallel for
        for (size_t i = 0; i < Size(); ++i)
            _M[i] += M[i];
#else
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::plus<T>());
#endif
#endif
        
        return *this;
        
    }

    
/**
 * @brief           ELementwise multiplication and assignment operator. i.e. this = m.
 *
 * @param  M        Added matrix.
 * @return          Result
 */
    template <class S> inline Matrix<T,P>&
    operator+=         (const Matrix<S,P>& M) {

        assert (_dim == M.Dim());

#ifdef EW_OMP
    #pragma omp parallel for
		for (size_t i = 0; i < Size(); ++i)
			_M[i] *= M[i];
#else
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::plus<T>());
#endif

		return *this;

    }


    /**
     * @brief           ELementwise addition with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
     * @return          Result
     */
    template <class S > inline Matrix<T,P>&
    operator+=          (const S s) {

    	T t = T (s);

#ifdef EW_OMP
    #pragma omp parallel for
		for (size_t i = 0; i < Size(); ++i)
			_M[i] += t;
#else
        std::transform (_M.begin(), _M.end(), _M.begin(), std::bind2nd(std::plus<T>(),t));
#endif
        return *this;

    }

    
    
    /**
     * @brief           Elementwise substruction of two matrices
     *
     * @param  M        Matrix substruent.
     */
    inline Matrix<T,P>
    operator-           (const Matrix<T,P>& M) const {

		Matrix<T,P> res = *this;
		return res -= M;
        
    }

    
    /**
     * @brief           Elementwise substruction of two matrices
     *
     * @param  M        Matrix substruent.
     */
    template <class S> inline Matrix<T,P>
    operator-           (const Matrix<S,P>& M) const {

		Matrix<T,P> res = *this;
		return res -= M;

    }


    /**
     * @brief           Elementwise subtraction all elements by a scalar
     *
     * @param  s        Scalar substruent.
     */
    template <class S>
    inline Matrix<T,P>
    operator-           (const S s) const {

		Matrix<T,P> res = *this;
		return res -= s;

    }

    
    
    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator-=         (const Matrix<T,P>& M) {

        assert (_dim == M._dim);

#if defined USE_VALARRAY
        _M -= M.Container();
#elif defined EXPLICIT_SIMD
        if (fp_type(_M[0]))
        	SSE::binary<T>(_M, M.Container(), SSE::sub<T>(), _M);
        else
        	for (size_t i = 0; i < Size(); ++i)
        		_M[i] -= M[i];
#else
#ifdef EW_OMP
    #pragma omp parallel for
        for (size_t i = 0; i < Size(); ++i)
            _M[i] -= M[i];
#else
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::minus<T>());
#endif
#endif

        return *this;

    }


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = m.
     *
     * @param  M        Added matrix.
     * @return          Result
     */
    template <class S> inline Matrix<T,P>&
    operator-=          (const Matrix<S,P>& M) {

        assert (_dim == M.Dim());

#ifdef EW_OMP
    #pragma omp parallel for
		for (size_t i = 0; i < Size(); ++i)
			_M[i] -= M[i];
#else
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::minus<T>());
#endif

        return *this;

    }

    /**
     * @brief           ELementwise substration with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Added scalar.
     * @return          Result
     */
    template <class S >
    inline Matrix<T,P>&
    operator-=          (const S s) {

		T t = T (s);

#ifdef EW_OMP
    #pragma omp parallel for
		for (size_t i = 0; i < Size(); ++i)
			_M[i] -= t;
#else
        std::transform (_M.begin(), _M.end(), _M.begin(), std::bind2nd(std::minus<T>(),t));
#endif

		return *this;

    }

    /**
     * @brief           Elementwise multiplication. i.e. this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>
    operator*          (const Matrix<T,P> &M) const {

		Matrix<T,P> res = *this;
		return res *= M;

    }


    /**
     * @brief           Elementwise multiplication. i.e. this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>
    operator*          (const Matrix<S,P> &M) const {

		Matrix<T,P> res = *this;
		return res *= M;

    }


    /**
     * @brief           Elementwise multiplication with a scalar. i.e. this * m.
     *
     * @param  s        Factor scalar
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>
    operator*          (const S s) const  {

        Matrix<T,P> res = *this;
        return res *= s;

    }



    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator*=         (const Matrix<T,P>& M) {

        assert (_dim == M._dim);

#if defined USE_VALARRAY
        _M *= M.Container();
#elif defined EXPLICIT_SIMD
        if (fp_type(_M[0]))
        	SSE::binary<T>(_M, M.Container(), SSE::mul<T>(), _M);
        else
        	for (size_t i = 0; i < Size(); ++i)
        		_M[i] *= M[i];
#else
#ifdef EW_OMP
    #pragma omp parallel for
        for (size_t i = 0; i < Size(); ++i)
            _M[i] *= M[i];
#else
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::multiplies<T>());
#endif
#endif

        return *this;

    }


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    template <class S> inline Matrix<T,P>&
    operator*=         (const Matrix<S,P>& M) {

        assert (_dim == M.Dim());

#ifdef EW_OMP
    #pragma omp parallel for
		for (size_t i = 0; i < Size(); ++i)
			_M[i] *= M[i];
#else
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::multiplies<T>());
#endif

		return *this;

    }


    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this *= s.
     *
     * @param  s        Factor scalar.
     * @return          Result
     */
    template <class S> inline Matrix<T,P>&
    operator*=         (const S s) {

#ifdef EW_OMP
    #pragma omp parallel for
		for (size_t i = 0; i < Size(); ++i)
			_M[i] += s;
#else
        std::transform (_M.begin(), _M.end(), _M.begin(), std::bind2nd(std::multiplies<T>(),s));
#endif
		return *this;

    }


    /**
     * @brief           Elementwise substruction of two matrices
     *
     * @param  M        Matrix substruent.
     */
    inline Matrix<T,P>
    operator/           (const Matrix<T,P>& M) const {

		Matrix<T,P> res = *this;
		return res /= M;

    }


    /**
     * @brief           Elelemtwise division by M.
     *
     * @param  M        The divisor.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>
    operator/          (const Matrix<S,P>& M) const {

		Matrix<T,P> res = *this;
		return res /= M;

    }


    /**
     * @brief           Elementwise division by scalar. i.e. this * m.
     *
     * @param  s        The divisor.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>
    operator/           (const S s) const {

		Matrix<T,P> res = *this;
		return res /= s;

	}


    /**
     * @brief           ELementwise multiplication and assignment operator. i.e. this = this .* M.
     *
     * @param  M        Factor matrix.
     * @return          Result
     */
    inline Matrix<T,P>&
    operator/=         (const Matrix<T,P>& M) {

        assert (_dim == M._dim);

#if defined USE_VALARRAY
        _M /= M.Container();
#elif defined EXPLICIT_SIMD
        if (fp_type(_M[0]))
        	SSE::binary<T>(_M, M.Container(), SSE::div<T>(), _M);
        else
        	for (size_t i = 0; i < Size(); ++i)
        		_M[i] /= M[i];
#else
#ifdef EW_OMP
    #pragma omp parallel for
        for (size_t i = 0; i < Size(); ++i)
            _M[i] /= M[i];
#else
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::divides<T>());
        #endif
#endif

        return *this;

    }


    /**
     * @brief           ELementwise division and assignment operator. i.e. this = this ./ M.
     *
     * @param  M        Divisor matrix.
     * @return          Result
     */
    template <class S> inline Matrix<T,P>&
    operator/=         (const Matrix<S,P> &M) {

        assert (_dim == M.Dim());

#ifdef EW_OMP
    #pragma omp parallel for
		for (size_t i = 0; i < Size(); ++i)
			_M[i] /= M[i];
#else
        std::transform (_M.begin(), _M.end(), M.Begin(), _M.begin(), std::divides<T>());
#endif
        return *this;

    }



    /**
     * @brief           ELementwise multiplication with scalar and assignment operator. i.e. this = m.
     *
     * @param  s        Divisor scalar.
     * @return          Result
     */
    template <class S>
    inline Matrix<T,P>&
    operator/=         (const S s) {

#ifdef EW_OMP
    #pragma omp parallel for
		for (size_t i = 0; i < Size(); ++i)
			_M[i] /= T(s);
#else
        std::transform (_M.begin(), _M.end(), _M.begin(), std::bind2nd(std::divides<T>(),s));
#endif

        return *this;

    }


    //@}


    /**
     * @name            Friend operators
     *                  Who doesn't need friends
     */

    //@{


    //--
    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator*  (const double s, const Matrix<T,P>& m) {
        return   m * s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator*  (const float s, const Matrix<T,P> &m) {
        return   m * s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator*  (const short s, const Matrix<T,P> &m) {
        return   m * s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator*  (const long s, const Matrix<T,P> &m) {
        return   m * s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator*  (const cxfl s, const Matrix<T,P> &m) {
        return   m * s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator*  (const cxdb s, const Matrix<T,P> &m) {
        return   m * s;
    }


    //--
    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator+  (const double s, const Matrix<T,P> &m) {
        return   m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator+  (const float s, const Matrix<T,P> &m) {
        return   m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator+  (const short s, const Matrix<T,P> &m) {
        return   m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator+  (const long s, const Matrix<T,P> &m) {
        return   m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator+  (const cxfl s, const Matrix<T,P> &m) {
        return   m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator+  (const cxdb s, const Matrix<T,P> &m) {
        return   m + s;
    }


    //--
    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator-  (const double s, const Matrix<T,P> &m) {
        return -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator-  (const float s, const Matrix<T,P> &m) {
        return -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator-  (const short s, const Matrix<T,P> &m) {
        return -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator-  (const long s, const Matrix<T,P> &m) {
        return   -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator-  (const cxfl s, const Matrix<T,P> &m) {
        return   -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator-  (const cxdb s, const Matrix<T,P> &m) {
        return   -m + s;
    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator/  (const double s, const Matrix<T,P> &m) {

        Matrix<T,P> res = m;
#ifdef USE_VALARRAY
		res.Container() = s / res.Container();
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < m.Size(); ++i)
            res[i] = s / res[i];
#endif
        return res;

    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator/  (const float s, const Matrix<T,P> &m) {

		Matrix<T,P> res = m;
#ifdef USE_VALARRAY
		res.Container() = s / res.Container();
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < m.Size(); ++i)
            res[i] = s / res[i];
#endif
		return res;

    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator/  (const cxfl s, const Matrix<T,P> &m) {

        Matrix<T,P> res = m;
#ifdef USE_VALARRAY
		res.Container() = s / res.Container();
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < m.Size(); ++i)
            res[i] = s / res[i];
#endif
        return res;

    }


    /**
     * @brief           Elementwise multiplication with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m * s
     */
    inline friend Matrix<T,P>
    operator/  (const cxdb s, const Matrix<T,P> &m) {

		Matrix<T,P> res = m;

#ifdef USE_VALARRAY
		res.Container() = s / res.Container();
#else
#ifdef EW_OMP
    #pragma omp parallel for
#endif
        for (size_t i = 0; i < m.Size(); ++i)
            res[i] = s / res[i];
#endif

		return res;

    }


    /**
     * @brief           Elementwise equality with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m == s
     */
    inline friend Matrix<cbool>
    operator== (const T s, const Matrix<T,P>& m) {
        return   m == s;
    }


    /**
     * @brief           Elementwise >= with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          m <= t
     */
    inline friend Matrix<cbool>
    operator>= (const T s, const Matrix<T,P>& m) {
        return   m <= s;
    }


    /**
     * @brief           Elementwise <= with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T<=M
     */
    inline friend Matrix<cbool>
    operator<= (const T s, const Matrix<T,P>& m) {
        return   m >= s;
    }


    /**
     * @brief           Elementwise unequality with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T!=M
     */
    inline friend Matrix<cbool>
    operator!= (const T s, const Matrix<T,P>& m) {
        return   m != s;
    }


    /**
     * @brief           Elementwise equality with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T+M
     */
    inline friend Matrix<cbool>
    operator>  (const T s, const Matrix<T,P>& m) {
        return   m <  s;
    }


    /**
     * @brief           Elementwise < with scalar (lhs)
     *
     * @param  s        Scalar lhs
     * @param  m        Matrix rhs
     * @return          T+M
     */
    inline friend Matrix<cbool>
    operator<  (const T s, const Matrix<T,P>& m) {
        return   m >  s;
    }


    /**
     * @brief           Elementwise equality with scalar (lhs)
     *
     * @param  mb       Scalar lhs
     * @param  m        Matrix rhs
     * @return          T+M
     */
    inline friend Matrix<T,P>
    operator&  (const Matrix<cbool>& mb, const Matrix<T,P>& m) {
        return   m & mb;
    }

    //@}


    /**
     * @brief           Elementwise equality, result[i] = (this[i] == m[i]). i.e. this == m
     *
     * @param  M        Comparing matrix.
     * @return          Hit list
     */
    inline Matrix<cbool>
    operator==          (const Matrix<T,P>& M) const {

        assert (Size() == M.Size());

        Matrix<cbool> res(_dim);
#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] = (_M[i] == M[i]) ? 1 : 0;

        return res;

    }


    /**
	 * @brief           Elementwise equality, result[i] = (this[i] == m[i]). i.e. this == m
	 *
	 * @param  M        Comparing matrix.
	 * @return          Hit list
	 */
    template<class S>
	inline Matrix<cbool>
	operator==          (const Matrix<S,P>& M) const {

        assert (Size() == M.Size());

		Matrix<cbool> res (_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] = (_M[i] == (T)M[i]) ? 1 : 0;

		return res;

     }

    /**
     * @brief           Scalar equality. result[i] = (this[i] == m).
     *
     * @param  s        Comparing scalar.
     * @return          Matrix of true where elements are equal s and false else.
     */
    inline Matrix<cbool>
    operator==          (const T s) const {

    	T t = (T) s;

        Matrix<cbool> res (_dim);

#ifdef EW_OMP
    #pragma omp parallel for
#endif
		for (size_t i = 0; i < Size(); ++i)
			res[i] =  (_M[i] == s) ? 1 : 0;

        return res;

    }
