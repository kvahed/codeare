#ifndef __CODEARE_FT_H__
#define __CODEARE_FT_H__

template<class T> class FT {
public:
	FT ();
	FT (const Params&);
	virtual ~FT ();
	virtual Matrix< std::complex<T> > Trafo (const Matrix< std::complex<T> >&) const;
	virtual Matrix< std::complex<T> > Adjoint (const Matrix< std::complex<T> >&) const;
	virtual Matrix< std::complex<T> > operator* (const Matrix< std::complex<T> >&) const;
	virtual Matrix< std::complex<T> > operator->* (const Matrix< std::complex<T> >&) const;
};

template<class T> inline Matrix<T> fftshift (const Matrix<T>&);
template<class T> inline Matrix< std::complex<T> > hannwindow (const Matrix<size_t>&, const T&);

template<class T> class DFT : public FT<T> {
public:
    DFT ();
	DFT (const Matrix<size_t>&);
	DFT (const Params&);
	virtual ~DFT ();
	virtual Matrix< std::complex<T> > Trafo (const Matrix< std::complex<T> >&) const;
	virtual Matrix< std::complex<T> > Adjoint (const Matrix< std::complex<T> >& m) const;
	virtual Matrix< std::complex<T> > operator* (const Matrix< std::complex<T> >&) const;
	virtual Matrix< std::complex<T> > operator->* (const Matrix< std::complex<T> >&) const;
};

template<class T> class NCSENSE : public FT<T> {
public:
	NCSENSE ();
	NCSENSE (const Params&);
	virtual ~NCSENSE ();
	void KSpace (const Matrix<T>&);
	void Weights (const Matrix<T>&);
	virtual Matrix< std::complex<T> > Trafo (const Matrix< std::complex<T> >&);
	virtual Matrix< std::complex<T> > Trafo (const Matrix< std::complex<T> >&, const Matrix< std::complex<T> >&, const bool&) const;
	virtual Matrix< std::complex<T> > Adjoint (const Matrix< std::complex<T> >&) const;
	virtual Matrix< std::complex<T> > Adjoint (const Matrix< std::complex<T> >&, const Matrix< std::complex<T> >&, const bool&) const;
	virtual Matrix< std::complex<T> > operator* (const Matrix< std::complex<T> >&) const;
	virtual Matrix< std::complex<T> > operator->* (const Matrix< std::complex<T> >&) const;
};

#endif
