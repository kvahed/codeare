/*
 * RangeParser.hpp
 *
 *  Created on: Mar 18, 2015
 *      Author: kvahed
 */

#ifndef RANGEPARSER_HPP_
#define RANGEPARSER_HPP_

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <ctype.h>

static inline std::vector<std::string>
Parse     (const std::string& str, const std::string& dlm) {

	assert (dlm.size() > 0);
	std::vector<std::string> sv;
	size_t  start = 0, end = 0;

	while (end != std::string::npos) {
		end = str.find (dlm, start);
		sv.push_back(str.substr(start, (end == std::string::npos) ?
				std::string::npos : end - start));
		start = ((end > (std::string::npos - dlm.size())) ?
				std::string::npos : end + dlm.size());
	}

	return sv;

}

enum RangeParseException {
	IMPROPER_RANGE_DECLARATION = 301,
	RANGE_DOES_NOT_FIT_MATRIX_DIMS,
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
 * @brief
 */
inline static Vector<size_t> Range () {
	return Vector<size_t>();
}
inline static Vector<size_t> Range (int i) {
	return Vector<size_t>(1,i);
}
inline static Vector<size_t> Range (int begin, int stride, int end) {
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
	if (stride>0)
		end++;
	else
		end--;
	Vector<size_t> view (std::ceil((float)(end-begin)/(float)stride));
	for (size_t i = 0; i < view.size(); ++i) {
		view[i] = begin;
		begin  += stride;
	}
	return view;
}
inline static Vector<size_t> Range (int begin, int end) {
	return Range(begin, 1, end);
}

#include <boost/tuple/tuple.hpp>



//template<class T> class Matrix;
#include "Vector.hpp"

inline static Vector<Vector<size_t> >
RangeParser (const std::string& range, const Vector<size_t>& dims) {
    std::string m_range (range);
    // remove white spaces
    m_range.erase(std::remove_if(m_range.begin(), m_range.end(), ::isspace), m_range.end());
    std::vector<std::string> ranges = Parse(m_range, ","); // Split ranges with commas.
    Vector<Vector<size_t> > view;
    
    for (size_t j = 0; j < ranges.size(); ++j) {
        std::vector<std::string> parts = Parse (ranges[j], ":");
        Vector<int> lims;
        
        if (ranges[j] == std::string(":")) {
            lims.push_back(0);
            lims.push_back(((ranges.size() == 1)?prod(dims):dims[j])-1);
        } else {
            for (size_t i = 0; i < parts.size(); ++i) {
                try {
                    if (parts[i].find(std::string("end")) != std::string::npos) {
                        lims.push_back(((ranges.size() == 1) ? prod(dims) : dims[j])-1);
                        parts[i] = parts[i].substr(3,parts[i].size()-3);
                        if (parts[i].size())
                            lims[i]+=boost::lexical_cast<int>(parts[i]);
                    }
                    else
                        lims.push_back(boost::lexical_cast<int>(parts[i]));
                } catch(const boost::bad_lexical_cast &) {
                    printf ("Failed to parse range, %s\n", ranges[j].c_str());
                    throw FAILED_TO_PARSE_RANGE;
                }
            }
            if (lims.size()==1)
            	lims.PushBack(lims[0]);
        }

        if (ranges.size() == 1) {
            view.PushBack(Range(lims[0],lims[1]));
            return view;
        }
        
        if (ranges.size() != dims.size()) {
            std::cout << "  ** ERROR - RangeParser **: "
            		"Range must cover equal number of dimensions " << std::endl;
            std::cout << "                            " << range <<
            		" does not match dimensions " << dims << std::endl;
            throw RANGE_DOES_NOT_FIT_MATRIX_DIMS;
        }
        
        if (parts.size() < 1 || parts.size() > 3) {
            printf ("Improper range declaration, \"%s\"\n", ranges[j].c_str());
            throw RANGE_MUST_CONSIST_OF_ONE_THROUGH_THREE_PARTS;
        } else if (parts.size() < 3) { // "x:y"
            if (lims[0]>dims[j] || lims[1]>dims[j]) {
                printf ("Exceeding matrix dimensions, \"%s\"\n", ranges[j].c_str());
                throw EXCEEDING_MATRIX_DIMENSIONS;
            }
            if (lims[0]>lims[1] || lims[0]<0 || lims[1]<0) {
                printf ("Improper range declaration, \"%s\"\n", ranges[j].c_str());
                throw EXCEEDING_MATRIX_DIMENSIONS;
            } else {
                view.PushBack(Range(lims[0],lims[1]));
            }
        } else if (parts.size() == 3) { // "x:y:z"
            throw YET_TO_BE_CODED;
        }
    }
    return view;
}

#endif /* RANGEPARSER_HPP_ */
