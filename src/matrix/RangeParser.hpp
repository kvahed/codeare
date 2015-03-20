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
		sv.push_back(str.substr(start, (end == std::string::npos) ? std::string::npos : end - start));
		start = ((end > (std::string::npos - dlm.size())) ? std::string::npos : end + dlm.size());
	}

	return sv;

}

enum RangeParseException {
	RANGE_DOES_NOT_FIT_MATRIX_DIMS
};

//template<class T> class Matrix;
#include "Vector.hpp"

inline static Vector<Vector<size_t> > RangeParser (const std::string& range, const Vector<size_t>& dims) {
    std::string m_range (range);
    m_range.erase(std::remove_if(m_range.begin(), m_range.end(), ::isspace), m_range.end()); // remove white spaces
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
                    printf ("Improper range declaration, %s\n", ranges[j].c_str());
                    assert (false);
                }
            }
        }

        if (ranges.size() == 1) {
            lims[1]++;
            view.PushBack(Vector<size_t>(lims[1]-lims[0]));
            for (size_t i = 0; i < view[0].size(); ++i)
                view[0][i] = lims[0] + i;
            return view;
        }
        
        if (ranges.size() != dims.size()) {
            std::cout << "  ** ERROR - RangeParser **: Range must cover equal number of dimensions " << std::endl;
            std::cout << "                            " << range << " does not match dimensions " << dims << std::endl;
            throw RANGE_DOES_NOT_FIT_MATRIX_DIMS;
        }
        
        if (parts.size() < 2 || parts.size() > 3) {
            printf ("Improper range declaration, %s\n", ranges[j].c_str());
            assert (false);
        } else if (parts.size() == 2) { // "x:y"
            lims[1]++;
            if (lims[0]>dims[j] || lims[1]>dims[j]) {
                std::cout << lims << std::endl;
                printf ("Exceeding range, \"%s\"\n", ranges[j].c_str());
                assert (false);
            }
            if (lims[0]>=lims[1] || lims[0]<0 || lims[1]<0) {
                printf ("Improper range declaration, %s\n", ranges[j].c_str());
                assert (false);
            } else {
                view.PushBack(Vector<size_t>(lims[1]-lims[0]));
                for (size_t i = 0; i < view[j].size(); ++i)
                    view[j][i] = lims[0] + i;
            }
        } else if (parts.size() == 3) { // "x:y:z"
            assert (false);
        }
    }
    return view;
}

#endif /* RANGEPARSER_HPP_ */
