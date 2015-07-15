#ifndef __CONST_NO_CONST_VIEW_HPP__
#define __CONST_NO_CONST_VIEW_HPP__

#include <type_traits>
#include <string>
#include <vector>

enum RangeParseException {
	IMPROPER_RANGE_DECLARATION = 301,
	RANGE_DOES_NOT_FIT_MATRIX_DIMS,
    CONCAT_MUST_CONSIST_OF_MINIMUM_ONE_PARTS,
	RANGE_MUST_CONSIST_OF_ONE_THROUGH_THREE_PARTS,
	FAILED_TO_PARSE_RANGE,
	EXCEEDING_MATRIX_DIMENSIONS,
	YET_TO_BE_CODED,
	NEGATIVE_BEGIN_INDEX,
	NEGATIVE_END_INDEX,
	NEGATIVE_STRIDE_REQUIRES_NEGATIV_RANGE,
	POSITIVE_STRIDE_REQUIRES_POSITIV_RANGE,
	STRIDE_MUST_NOT_BE_ZERO
};

template<bool is_const> class Range {

public:

	inline Range () {}

	inline Range (const size_t& begend) { HandleSingleInput(begend); }
    
	inline Range (const size_t& begin, const size_t& end) { HandleTwoInputs(begin, end); }
    
	inline Range (const size_t& begin, const size_t& stride, const size_t& end) {
		HandleThreeInputs(begin, stride, end);
	}
    
	inline Range (const Vector<size_t>& v) { _idx = v; }

	inline Range (const std::string& rs) { ParseRange(rs); }

	inline Range (const char* rcs) { ParseRange(std::string(rcs)); }

	inline void ParseRange (std::string rs) {
		// remove white spaces
		rs.erase(std::remove_if(rs.begin(), rs.end(), ::isspace), rs.end());
		// split delimited by kommas
		std::vector<std::string> concat = Parse(rs,",");
		if (concat.empty())
			throw CONCAT_MUST_CONSIST_OF_MINIMUM_ONE_PARTS;
		for (size_t i = 0; i < concat.size(); ++i) {
			std::vector<std::string> parts = Parse(concat[i],":");
			switch (parts.size()) {
			case 1:
				HandleSingleInput(boost::lexical_cast<int>(parts[0]));
				break;
			case 2:
				HandleTwoInputs(boost::lexical_cast<int>(parts[0]),
								boost::lexical_cast<int>(parts[1]));
				break;
			case 3:
				HandleThreeInputs(boost::lexical_cast<int>(parts[0]),
								  boost::lexical_cast<int>(parts[1]),
								  boost::lexical_cast<int>(parts[2]));
				break;
			default:
				throw RANGE_MUST_CONSIST_OF_ONE_THROUGH_THREE_PARTS;
				break;
			}
		}
	}
    
	virtual ~Range() {}

	inline void Reset (const int& begin) {
		_idx.Clear();
		HandleSingleInput (begin);
	}
    
	inline void Reset (const int& begin, const int& end) {
		_idx.Clear();
		HandleTwoInputs (begin, end);
	}
    
	inline void Reset (const int& begin, const int& stride, const int& end) {
		_idx.Clear();
		HandleThreeInputs (begin, stride, end);
	}
    
	inline size_t Size() const { return _idx.size(); }
    
	inline bool IsSingleton() const { return (_idx.size()==1);}
    
	inline size_t operator[] (const size_t& i) const { return _idx[i]; }
    
private:

	inline void HandleSingleInput (const int& pos) {
		if (pos < 0)
			throw  NEGATIVE_BEGIN_INDEX;
		_idx.push_back(pos);
	}

	inline void HandleThreeInputs (const int& begin, const int& stride, const int& end) {
		if (begin < 0) {
			printf ("NEGATIVE_BEGIN_INDEX\n");
			throw NEGATIVE_BEGIN_INDEX;
		}
		if (end < 0) {
			printf ("NEGATIVE_END_INDEX\n");
			throw NEGATIVE_END_INDEX;
		}
		if (stride == 0) {
			printf ("STRIDE_MUST_NOT_BE_ZERO\n");
			throw STRIDE_MUST_NOT_BE_ZERO;
			}
		if (stride > 0 && end < begin) {
			printf ("POSITIVE_STRIDE_REQUIRES_POSITIV_RANGE\n");
			throw POSITIVE_STRIDE_REQUIRES_POSITIV_RANGE;
		}
		if (stride < 0 && end > begin) {
			printf ("NEGATIVE_STRIDE_REQUIRES_NEGATIV_RANGE\n");
			throw NEGATIVE_STRIDE_REQUIRES_NEGATIV_RANGE;
		}
		size_t cur = _idx.size();
		_idx.resize(cur+std::floor(((float)end-(float)begin)/(float)stride)+1);
		for (size_t i = cur; i < _idx.size(); ++i)
			_idx[i] = begin + i*stride;
	}

	inline void HandleTwoInputs (const int& begin, const int& end) {
		if (begin < 0) {
			printf ("NEGATIVE_BEGIN_INDEX\n");
			throw NEGATIVE_BEGIN_INDEX;
		}
		if (end < 0) {
			printf ("%d: NEGATIVE_END_INDEX\n", end);
			throw NEGATIVE_END_INDEX;
		}
		if (end < begin) {
			printf ("POSITIVE_STRIDE_REQUIRES_POSITIV_RANGE\n");
			throw STRIDE_MUST_NOT_BE_ZERO;
		}
		size_t cur = _idx.size();
		_idx.resize(cur+end-begin+1);
		for (size_t i = cur; i < _idx.size(); ++i)
			_idx[i] = begin + i;
	}

	friend std::ostream& operator<< (std::ostream &os, const Range& r) {
		return os << r._idx;
	}
	Vector<size_t> _idx;
};

typedef Range<true> CR;
typedef Range<false> R;

template <class T, bool is_const> class ConstNoConstView;
template <class T, paradigm P=SHM> class Matrix;

template<class T, paradigm P=SHM>
class MatrixType {
public:
    typedef ConstNoConstView<T,true> ConstView;
    typedef ConstNoConstView<T,false> View;

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
	virtual View operator() (R r) { return View(); }
	virtual ConstView operator() (const CR r) const { return ConstView(); }
	virtual View operator() (const R r, const size_t& n) { return View(); }
	virtual ConstView operator() (const CR r, const size_t& n) const { return ConstView(); }
	virtual View operator() (R r0, R r1) { return View(); }
	virtual ConstView operator() (CR r0, CR r1) const { return ConstView(); }
	virtual View operator() (R r0, R r1, R r2) { return View(); }
	virtual ConstView operator() (CR r0, CR r1, CR r2) const { return ConstView(); }
	virtual View operator() (R r0, R r1, R r2, R r3) { return View(); }
	virtual ConstView operator() (CR r0, CR r1, CR r2, CR r3) const { return ConstView(); }
	virtual View operator() (R r0, R r1, R r2, R r3, R r4) { return View(); }
	virtual ConstView operator() (CR r0, CR r1, CR r2, CR r3, CR r4) const { return ConstView(); }
private:
    Vector<T> _M;
	Vector<size_t> _dim;
};

template<class T, bool is_const = true> class ConstNoConstView : public MatrixType<T> {
public:
    
    typedef typename std::conditional<is_const, const Matrix<T>, Matrix<T> >::type MatrixTypeType;
    typedef typename std::conditional<is_const, const T, T>::type Type;
    
    inline ConstNoConstView () : _matrix(0) {}
    inline ConstNoConstView (MatrixTypeType* matrix, Vector<Range<is_const> >& range) :
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

        while (_range.back().IsSingleton())
            _range.PopBack();

    }
    
    operator Matrix<T>() const {
        Matrix<T> res (_dim);
        for (size_t i = 0; i < Size(); ++i)
            res[i] = *(_pointers[i]);
            return res;
    }
    
    inline ConstNoConstView& operator= (const Matrix<T>& M) {
        assert (Size() == M.Size());
        for (size_t i = 0; i < Size(); ++i)
            *(_pointers[i]) = M[i];
        return *this;
    }
        
    inline virtual const T& operator[] (const size_t& pos) const {
        assert(pos < Size());
        return *(_pointers[pos]);
    }

    template<class S>
    inline Matrix<T> operator* (const MatrixType<S>& d) const {
        assert (Size() == d.Size());
        Matrix<T> M(_dim);
        for (size_t i = 0; i < Size(); ++i)
            M[i] = *(_pointers[i])*d[i];
        return M;
    }
    
    inline virtual Matrix<T> operator/ (const MatrixTypeType& d) const {
        Matrix<T> M(_dim);
        for (size_t i = 0; i < Size(); ++i)
            M[i] = *(_pointers[i])/d[i];
        return M;
    }
    
    inline virtual Matrix<T> operator+ (const MatrixTypeType& d) const {
        Matrix<T> M(_dim);
        for (size_t i = 0; i < Size(); ++i)
            M[i] = *(_pointers[i])+d[i];
        return M;
    }

    inline virtual Matrix<T> operator- (const MatrixTypeType& d) const {
        Matrix<T> M(_dim);
        for (size_t i = 0; i < Size(); ++i)
            M[i] = *(_pointers[i])-d[i];
        return M;
    }

    virtual ~ConstNoConstView () { _matrix = 0; }
    
    inline ConstNoConstView& operator= (const ConstNoConstView<T,true>& v) {
        Matrix<T>& lhs = *_matrix;
        const Matrix<T>& rhs = *(v._matrix);
        assert (_nsdims.size() == v._nsdims.size());
        for (size_t i = 0; i < _nsdims.size(); ++i)
            assert(_range[i].Size()==v._range[i].Size());
        for (size_t i = 0; i < _pointers.size(); ++i)
            *_pointers[i] = *(v._pointers)[i];
        return *this;
    }
    inline ConstNoConstView& operator= (const Type& t) {
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
    friend std::ostream& operator<< (std::ostream &os, const ConstNoConstView& r) {
        os << "(";
        for (size_t i = 0; i < r._range.size(); ++i) {
            if (i)
                os << ",";
            os << r._range[i];
        }
        return os << ")" << std::endl;
    }
};

#endif
