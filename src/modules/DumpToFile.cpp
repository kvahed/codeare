/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but 
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301  USA
 */

#include "DumpToFile.h"

using namespace RRStrategy;

DumpToFile::DumpToFile () {}

error_code 
DumpToFile::Process () {

	std::stringstream fname;
	const char* uid = Attribute ("UID");

	if (uid == 0  ||  uid == "")
		uid = "unspecified";

	fname << uid << "_raw.bin";
	m_raw.dump    (fname.str().c_str());

	fname.str("");

	fname << uid << "_helper.bin";
	//m_helper.dump (fname.str().c_str());

	fname.str("");

	fname << uid << "_pixel.bin";
	//m_pixel.dump (fname.str().c_str());

	fname.str("");

	fname << uid << "_config.xml";
	DumpConfig (fname.str().c_str());

	return RRSModule::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new DumpToFile;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

