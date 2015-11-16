/*
 *  This file is part of rawpp.
 *
 *  rawpp is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  rawpp is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with rawpp. If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 17 Jul, 2014
 *      Author: Kaveh Vahedipour
 */

#ifndef _SYNGO_FILE_HPP_
#define _SYNGO_FILE_HPP_

#include "SyngoHeader.hpp"
//#include "Timer.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#ifdef USE_IN_MATLAB
#include <mex.h>
#endif
#include <algorithm>

static const float kB = 1.f/1024.f, iMB = 1.f/1.31072f, iGB = 1.f/1.073741824f;
static uint32_t header_len, data_len;

typedef std::complex<float> raw;

/**
 * @brief Helper. Get maximum of every value in lhs.
 */
inline static std::vector<uint32_t>& _max (std::vector<uint32_t>& dims, const uint16_t slc[14]) {
    if (slc[ 0] > dims[ 2]) dims[ 2] = slc[ 0]; if (slc[ 1] > dims[ 3]) dims[ 3] = slc[ 1];
    if (slc[ 2] > dims[ 4]) dims[ 4] = slc[ 2]; if (slc[ 3] > dims[ 5]) dims[ 5] = slc[ 3];
    if (slc[ 4] > dims[ 6]) dims[ 6] = slc[ 4]; if (slc[ 5] > dims[ 7]) dims[ 7] = slc[ 5];
    if (slc[ 6] > dims[ 8]) dims[ 8] = slc[ 6]; if (slc[ 7] > dims[ 9]) dims[ 9] = slc[ 7];
    if (slc[ 8] > dims[10]) dims[10] = slc[ 8]; if (slc[ 9] > dims[11]) dims[11] = slc[ 9];
    if (slc[10] > dims[12]) dims[12] = slc[10]; if (slc[11] > dims[13]) dims[13] = slc[11];
    if (slc[12] > dims[14]) dims[14] = slc[12]; if (slc[13] > dims[15]) dims[15] = slc[13];
    return dims;
}

/**
 * @brief Helper. Increment all values but first 2 in v.
 */
inline static std::vector<uint32_t>& _raise_one (std::vector<uint32_t>& v) {
    if (v[0]==0)
        v[0] = 1;
    if (v[1]==0)
        v[1] = 1;
    if (v.size()>2)
        for (size_t i = 2; i < v.size(); ++i)
            v[i]++;
    return v;
}

/**
 * @brief allocate individual left hand sides
 */
//inline static

#include "IOFile.hpp"

namespace codeare {
namespace matrix  {
namespace io      {

/**
 * @brief Base class for Syngo MR file readers
 */
class SyngoFile : public IOFile {
    
public:
    
#ifdef USE_IN_MATLAB
    /**
     * @brief Construct with filename and left hand side 
     */
    SyngoFile (const std::string& fname, int nlhs, mxArray *lhs[]) :
        _fname(fname), _header_len(0), _nlhs(nlhs), _lhs(lhs), _status(0) {
        prtmsg ("   Opening %s ...\n", _fname.c_str());
        std::string line;
        _file.open (fname.c_str(), std::ios::in|std::ios::binary);
        if (!_file.is_open()) {
            prterr ("     FAILED! Unable to open file\n");
            _status = 1;
        }
        prtmsg("     done.\n");
    }
#else
    /**
     * @brief Construct with filename and left hand side
     */
    SyngoFile (const std::string& fname, const IOMode mode,
			  const Params& params, const bool verbosity) :
        _fname(fname), _header_len(0), _nlhs(0), _lhs(0), _status(0) {
        prtmsg ("   Opening %s ...\n", _fname.c_str());
        std::string line;
        _file.open (fname.c_str(), std::ios::in|std::ios::binary);
        if (!_file.is_open()) {
            prterr ("     FAILED! Unable to open file\n");
            _status = 1;
        }
        prtmsg("     done.\n");
    }
#endif

    /**
     * @brief Close file
     */
    virtual ~SyngoFile () {
        this->Close();
    }
    
    /**
     * @brief Get status 
     */
    int Status () const {
        return _status;
    }
    
    /**
     * @brief Close file 
     */
    void Close () {
        if (_file.is_open()) {
            prtmsg ("   Closing %s ... \n", _fname.c_str());
            _file.close();
            prtmsg ("     done.\n");
        }
    }
    
    /**
     * @brief Digest file
     */
    virtual void Digest () = 0;
    

protected:
    
    int _status ;          /**< Well being of the reader */
    
    std::ifstream _file;   /**< File */
    std::string _fname;    /**< File name */
    
    uint32_t _header_len;  /**< Header length */
    int _nlhs;             /**< # Left hand sides */
#ifdef USE_IN_MATLAB
    mxArray** _lhs;        /**< Left hand side */
#else
    void* _lhs;
#endif

};

}}}
#endif
