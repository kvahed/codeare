#ifndef __CODEARE_MATRIX_H__
#define __CODEARE_MATRIX_H__

template <class T, paradigm P = SHM>
class Matrix  : public SmartObject {
    
public:

    Matrix ();
    Matrix (const std::vector<size_t>& dim);
    Matrix (const std::vector<size_t>& dim);
    Matrix (const std::vector<size_t>& dim, const std::vector<float>& res);
    Matrix (const size_t dim[INVALID_DIM], const float res[INVALID_DIM]);
    Matrix (const size_t n);
    Matrix (const size_t m, const size_t n);
    Matrix (const size_t m, const size_t n, const size_t k);
    Matrix (const size_t col, const size_t lin, const size_t cha, const size_t set,
            const size_t eco, const size_t phs, const size_t rep, const size_t seg,
            const size_t par, const size_t slc, const size_t ida, const size_t idb,
            const size_t idc, const size_t idd, const size_t ide, const size_t ave);
    Matrix (const Matrix<T,P> &M);
	virtual ~Matrix();
    T  operator[] (const size_t& p) const;
    T& operator[] (const size_t& p);
    T  operator() (const size_t& p) const;
    T& operator() (const size_t& p);
    T  operator() (const size_t x, const size_t y) const;
    T& operator() (const size_t x, const size_t y);
    T  operator() (const size_t x, const size_t y, const size_t z) const;
    T& operator() (const size_t x, const size_t y, const size_t z);
    T  operator() (const size_t col, const size_t lin, const size_t cha, const size_t set,
                   const size_t eco, const size_t phs, const size_t rep, const size_t seg,
                   const size_t par, const size_t slc, const size_t ida, const size_t idb,
                   const size_t idc, const size_t idd, const size_t ide, const size_t ave) const;
    T& operator() (const size_t col, const size_t lin, const size_t cha, const size_t set,
                   const size_t eco, const size_t phs, const size_t rep, const size_t seg,
                   const size_t par, const size_t slc, const size_t ida, const size_t idb,
                   const size_t idc, const size_t idd, const size_t ide, const size_t ave);

};

#endif
