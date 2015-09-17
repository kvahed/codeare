#ifndef __VIEW_HPP__
#define __VIEW_HPP__

#include <type_traits>
#include <string>
#include <vector>

#include "Range.hpp"

template <class T, bool is_const> class View;
template <class T, paradigm P=SHM> class Matrix;

template<class T, paradigm P=SHM>
class MatrixType {
public:
    typedef View<T,true> RHSView;
    typedef View<T,false> LHSView;

	virtual ~MatrixType() {}
    virtual const T& operator[] (const size_t& n) const { return _M[0]; }
    virtual size_t Size () const { return 0; }
	virtual Matrix<T> operator* (const MatrixType<T>& d) const { return Matrix<T>(); }
	virtual Matrix<T> operator/ (const MatrixType<T>& d) const { return Matrix<T>(); }
    virtual Matrix<T> operator+ (const MatrixType<T>& d) const { return Matrix<T>(); }
	virtual Matrix<T> operator- (const MatrixType<T>& d) const { return Matrix<T>(); }
	virtual Matrix<T> operator= (const MatrixType<T>& d) const { return Matrix<T>(); }
    virtual size_t NDim() const  { return 0; }
	virtual const Vector<size_t>& Dim() const { return _dim;  }
	virtual size_t Dim(const size_t& d) const { return _dim[0]; }
	virtual LHSView operator() (R r) { return LHSView(); }
	virtual RHSView operator() (const CR r) const { return RHSView(); }
	virtual LHSView operator() (const R r, const size_t& n) { return LHSView(); }
	virtual RHSView operator() (const CR r, const size_t& n) const { return RHSView(); }
	virtual LHSView operator() (R r0, R r1) { return LHSView(); }
	virtual RHSView operator() (CR r0, CR r1) const { return RHSView(); }
	virtual LHSView operator() (R r0, R r1, R r2) { return LHSView(); }
	virtual RHSView operator() (CR r0, CR r1, CR r2) const { return RHSView(); }
	virtual LHSView operator() (R r0, R r1, R r2, R r3) { return LHSView(); }
	virtual RHSView operator() (CR r0, CR r1, CR r2, CR r3) const { return RHSView(); }
	virtual LHSView operator() (R r0, R r1, R r2, R r3, R r4) { return LHSView(); }
	virtual RHSView operator() (CR r0, CR r1, CR r2, CR r3, CR r4) const { return RHSView(); }
private:
    Vector<T> _M;
	Vector<size_t> _dim;
};

template<class T, bool is_const = true> class View : public MatrixType<T> {
public:
    
    typedef typename std::conditional<is_const, const Matrix<T>, Matrix<T> >::type MatrixTypeType;
    typedef typename std::conditional<is_const, const T, T>::type Type;
    
    inline View () : _matrix(0) {}

    inline View (MatrixTypeType* matrix, Vector<Range<is_const> >& range) :
        _matrix(matrix), _range(range) {
        assert (_range.size());
        for (size_t i = 0; i < _range.size(); ++i) {
            if (!_range[i].IsSingleton()) {
                if (_range[i].Size() == 0) 
                    _range[i].Reset (0,_matrix->Dim(i)-1);
                _nsdims.push_back(i);
            }
        }
        
        if (_range.size() == 1)
            for (size_t n0 = 0; n0 < _range[0].Size(); ++n0)
                _pointers.push_back(matrix->Ptr()+_range[0][n0]);
        else if (_range.size() == 2)
            for (size_t n1 = 0; n1 < _range[1].Size(); ++n1)
                for (size_t n0 = 0; n0 < _range[0].Size(); ++n0)
                    _pointers.push_back(&((*_matrix)(_range[0][n0],_range[1][n1])));
        else if (_range.size() == 3)
            for (size_t n2 = 0; n2 < _range[2].Size(); ++n2)
                for (size_t n1 = 0; n1 < _range[1].Size(); ++n1)
                    for (size_t n0 = 0; n0 < _range[0].Size(); ++n0)
                        _pointers.push_back(&((*_matrix)(_range[0][n0],_range[1][n1],_range[2][n2])));
        else if (_range.size() == 4)
            for (size_t n3 = 0; n3 < _range[3].Size(); ++n3)
                for (size_t n2 = 0; n2 < _range[2].Size(); ++n2)
                    for (size_t n1 = 0; n1 < _range[1].Size(); ++n1)
                        for (size_t n0 = 0; n0 < _range[0].Size(); ++n0)
                            _pointers.push_back(&((*_matrix)
                                                  (_range[0][n0],_range[1][n1],
                                                   _range[2][n2],_range[3][n3])));
        else if (_range.size() == 5)
            for (size_t n4 = 0; n4 < _range[4].Size(); ++n4)
                for (size_t n3 = 0; n3 < _range[3].Size(); ++n3)
                    for (size_t n2 = 0; n2 < _range[2].Size(); ++n2)
                        for (size_t n1 = 0; n1 < _range[1].Size(); ++n1)
                            for (size_t n0 = 0; n0 < _range[0].Size(); ++n0)
                                _pointers.push_back(&((*_matrix)
                                                      (_range[0][n0],_range[1][n1],
                                                       _range[2][n2],_range[3][n3],
                                                       _range[4][n4])));
        else if (_range.size() == 6)
            for (size_t n5 = 0; n5 < _range[5].Size(); ++n5)
                for (size_t n4 = 0; n4 < _range[4].Size(); ++n4)
                    for (size_t n3 = 0; n3 < _range[3].Size(); ++n3)
                        for (size_t n2 = 0; n2 < _range[2].Size(); ++n2)
                            for (size_t n1 = 0; n1 < _range[1].Size(); ++n1)
                                for (size_t n0 = 0; n0 < _range[0].Size(); ++n0)
                                    _pointers.push_back(&((*_matrix)
                                                          (_range[0][n0],_range[1][n1],
                                                           _range[2][n2],_range[3][n3],
                                                           _range[4][n4],_range[5][n5])));
        else if (_range.size() == 7)
            for (size_t n6 = 0; n6 < _range[6].Size(); ++n6)
                for (size_t n5 = 0; n5 < _range[5].Size(); ++n5)
                    for (size_t n4 = 0; n4 < _range[4].Size(); ++n4)
                        for (size_t n3 = 0; n3 < _range[3].Size(); ++n3)
                            for (size_t n2 = 0; n2 < _range[2].Size(); ++n2)
                                for (size_t n1 = 0; n1 < _range[1].Size(); ++n1)
                                    for (size_t n0 = 0; n0 < _range[0].Size(); ++n0)
                                        _pointers.push_back(&((*_matrix)
                                                              (_range[0][n0],_range[1][n1],
                                                               _range[2][n2],_range[3][n3],
                                                               _range[4][n4],_range[5][n5],
                                                               _range[6][n6])));
        else if (_range.size() == 8)
            for (size_t n7 = 0; n7 < _range[7].Size(); ++n7)
                for (size_t n6 = 0; n6 < _range[6].Size(); ++n6)
                    for (size_t n5 = 0; n5 < _range[5].Size(); ++n5)
                        for (size_t n4 = 0; n4 < _range[4].Size(); ++n4)
                            for (size_t n3 = 0; n3 < _range[3].Size(); ++n3)
                                for (size_t n2 = 0; n2 < _range[2].Size(); ++n2)
                                    for (size_t n1 = 0; n1 < _range[1].Size(); ++n1)
                                        for (size_t n0 = 0; n0 < _range[0].Size(); ++n0)
                                            _pointers.push_back(&((*_matrix)
                                                                  (_range[0][n0],_range[1][n1],
                                                                   _range[2][n2],_range[3][n3],
                                                                   _range[4][n4],_range[5][n5],
                                                                   _range[6][n6],_range[7][n7])));
        
        for (auto it = _range.begin(); it != _range.end();) {
            _dim.push_back(it->Size());
            if(it->IsSingleton())
                it = _range.Erase(it);
            else
                ++it;
        }

    }
    
    operator Matrix<T>() const {
        Matrix<T> res (_dim);
        for (size_t i = 0; i < Size(); ++i)
            res[i] = *(_pointers[i]);
            return res;
    }
    
    inline View& operator= (const Matrix<T>& M) {
        assert (Size() == M.Size());
        for (size_t i = 0; i < Size(); ++i)
            *(_pointers[i]) = M[i];
        return *this;
    }
        
    inline virtual const T& operator[] (const size_t& pos) const {
        assert(pos < Size());
        return *(_pointers[pos]);
    }

    template<class S> inline Matrix<T> operator* (const MatrixType<S>& d) const {
        assert (Size() == d.Size());
        Matrix<T> M(_dim);
        for (size_t i = 0; i < Size(); ++i)
            M[i] = *(_pointers[i])*d[i];
        return M;
    }
    template<class S> inline MatrixType<T>& operator*= (const MatrixType<S>& d)  {
        assert (Size() == d.Size());
        for (size_t i = 0; i < Size(); ++i)
            *(_pointers[i]) *= d[i];
        return *this;
    }
    
    inline virtual Matrix<T> operator/ (const MatrixTypeType& d) const {
        Matrix<T> M(_dim);
        for (size_t i = 0; i < Size(); ++i)
            M[i] = *(_pointers[i])/d[i];
        return M;
    }
    template<class S> inline MatrixType<T>& operator/= (const MatrixType<S>& d)  {
        assert (Size() == d.Size());
        for (size_t i = 0; i < Size(); ++i)
            *(_pointers[i]) /= d[i];
        return *this;
    }
    
    inline virtual Matrix<T> operator+ (const MatrixTypeType& d) const {
        Matrix<T> M(_dim);
        for (size_t i = 0; i < Size(); ++i)
            M[i] = *(_pointers[i])+d[i];
        return M;
    }
    template<class S> inline MatrixType<T>& operator+= (const MatrixType<S>& d)  {
        assert (Size() == d.Size());
        for (size_t i = 0; i < Size(); ++i)
            *(_pointers[i]) += d[i];
        return *this;
    }

    inline virtual Matrix<T> operator- (const MatrixTypeType& d) const {
        Matrix<T> M(_dim);
        for (size_t i = 0; i < Size(); ++i)
            M[i] = *(_pointers[i])-d[i];
        return M;
    }
    template<class S> inline MatrixType<T>& operator-= (const MatrixType<S>& d)  {
        assert (Size() == d.Size());
        for (size_t i = 0; i < Size(); ++i)
            *(_pointers[i]) -= d[i];
        return *this;
    }

    virtual ~View () { _matrix = 0; }
    
    inline View& operator= (const View<T,true>& v) {
        Matrix<T>& lhs = *_matrix;
        const Matrix<T>& rhs = *(v._matrix);
        assert (_nsdims.size() == v._nsdims.size());
        for (size_t i = 0; i < _nsdims.size(); ++i)
            assert(_range[i].Size()==v._range[i].Size());
        for (size_t i = 0; i < _pointers.size(); ++i)
            *_pointers[i] = *(v._pointers)[i];
        return *this;
    }
    inline View& operator= (const Type& t) {
        assert (_matrix);
        for (size_t i = 0; i < _pointers.size(); ++i)
            *_pointers[i] = t;
        return *this;
    }
    
    inline Range<is_const>& Rng() { return _range; }
    inline virtual size_t Size() const {return _pointers.size();}
    inline virtual size_t Dim (const size_t& i) const { assert (i<_dim.size()); return _dim[i];}
    inline virtual const Vector<size_t>& Dim() const { return _dim; }
    virtual size_t NDim() const  { return _dim.size(); }
    
    MatrixTypeType* _matrix;
    Vector<Range<is_const> > _range;
    Vector<Type*> _pointers;
    Vector<size_t> _nsdims;
    Vector<size_t> _dim;
    
private:
    friend std::ostream& operator<< (std::ostream &os, const View& r) {
        os << "(";
        for (size_t i = 0; i < r._range.size(); ++i) {
            if (i)
                os << ",";
            os << r._range[i];
        }
        return os << ")" << std::endl;
    }
};

#endif // __VIEW_HPP__
