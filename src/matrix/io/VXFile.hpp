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
 *  Created on: Jul, 2014
 *      Author: Kaveh Vahedipour
 */

#ifndef _VX_FILE_HPP_
#define _VX_FILE_HPP_

#include "Matrix.hpp"
#include "IOFile.hpp"
#include "VBFile.hpp"
#include "VDFile.hpp"

enum idea_version {IDEA_VB, IDEA_VD};

/**
 * @brief Dispatcher for VB/VD binary versions
 */
namespace codeare {
namespace matrix {
namespace io {


class VXFile : public IOFile {

public:

    /**
     * @brief Construct with file name and output structures
     */
#ifdef USE_MATLAB
    VXFile (const std::string& fname, int nlhs = 0, mxArray *lhs[] = 0) {
        idea_version iv = CheckVersion (fname);
        _context = (iv == IDEA_VD) ?
            (SyngoFile*) new VD::VDFile(fname, nlhs, lhs) :
            (SyngoFile*) new VB::VBFile(fname, nlhs, lhs);
    }
#else
    /**
     * @brief Construct with filename and left hand side
     */
    VXFile (const std::string& fname, const IOMode mode = READ,
			  const Params& params = Params(), const bool verbosity = false) : _context(0) {
        idea_version iv = CheckVersion (fname);
        if (iv == IDEA_VD)
        	_context = (SyngoFile*) new VD::VDFile(fname, mode, params, verbosity);
        else
        	_context = (SyngoFile*) new VB::VBFile(fname, mode, params, verbosity);
    }
#endif

    /**
     * @brief Cleanup memory and close file
     */
    virtual ~VXFile () {
        delete (_context);
    }

    /**
     * @return Current status
     */
    int Status () {
        return _context->Status();
    }

    /**
     * @brief Digest ingredients
     */
    void Digest () {
        _context->Digest();
    }

    /**
     * @brief Close file
     */
    void Close () {
        _context->Close();
    }
    
    template<class T> Matrix<T> Read (const std::string& data_name = "meas") {
    	if (_version == IDEA_VB)
    		return ((VB::VBFile*)_context)->Read<T>(data_name);
    	else
    		return ((VD::VDFile*)_context)->Read<T>(data_name);
    }

    template<class T> Matrix<T> Read  (const TiXmlElement *) {
        throw 0;
    }
        
private:

    idea_version CheckVersion (const std::string& fname) {
        std::ifstream file (fname.c_str());
        uint32_t x[2];
        file.read ((char*)x, 2*sizeof(uint32_t));
        file.close();
        return (x[0] == 0 && x[1] <= 64) ? IDEA_VD : IDEA_VB;
    }
    
    /**
     * @brief Version context
     */
    SyngoFile* _context;
    idea_version _version;
};

}}}

#endif
