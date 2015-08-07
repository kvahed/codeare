#ifndef __RANGE_HPP__
#define __RANGE_HPP__

#include "Vector.hpp"

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

/**
 * @brief Index range
 */
template<bool is_const> class Range {

public:

    /**
     * @brief default constructor
     */
	inline Range () {}

    /**
     * @brief construct as single index
     *        i.e. A(:,6) is A(R(),R(5))
     */
	inline Range (const size_t& begend) { HandleSingleInput(begend); }

    /**
     * @brief construct with begin and end
     *        i.e. A(6:8,:) is A(R(5:7),R())
     */
	inline Range (const size_t& begin, const size_t& end) { HandleTwoInputs(begin, end); }

    /**
     * @brief construct with begin, stride and end
     *        i.e. A(:,6:-1:4) translates to A(R(),R(5,-1,3))
     */
    inline Range (const size_t& begin, const size_t& stride, const size_t& end) {
		HandleThreeInputs(begin, stride, end);
	}

    /**
     * @brief construct with index vector
     *        (i.e. A(:,v) translates to A(R(),R(v)))
     */
	inline Range (const Vector<size_t>& v) { _idx = v; }

    /**
     * @brief construct with matlab-like string
     *         (i.e. A(:,6:-1:4) translates to A(":","5:-1:3")
     */
	inline Range (const std::string& rs) { ParseRange(rs); }

    /**
     * @brief construct with matlab-like string
     * @see Range (const std::string& rs)
     */
	inline Range (const char* rcs) { ParseRange(std::string(rcs)); }
    
    /**
     * @brief default destructor
     */
    virtual ~Range() {}

    /**
     * @brief reset range
     * @see Range(const int& begin)
     */
	inline void Reset (const int& begin) {
		_idx.Clear();
		HandleSingleInput (begin);
	}

    /**
     * @brief reset range
     * @see Range(const int& begin, const int& end)
     */
    inline void Reset (const int& begin, const int& end) {
		_idx.Clear();
		HandleTwoInputs (begin, end);
	}

    /**
     * @brief reset range
     * @see Range(const int& begin, const int& stride, const int& end)
     */
    inline void Reset (const int& begin, const int& stride, const int& end) {
		_idx.Clear();
		HandleThreeInputs (begin, stride, end);
	}

    /**
     * @brief Get size
     * @return Size
     */
	inline size_t Size() const { return _idx.size(); }

    /**
     * @brief Is this range a singleton
     * @return Singleton or not
     */
	inline bool IsSingleton() const { return (_idx.size()==1);}

    /**
     * @brief Get i-th index in range
     * @return I-th index in range
     */
	inline size_t operator[] (const size_t& i) const { return _idx[i]; }
    
private:

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

#endif
