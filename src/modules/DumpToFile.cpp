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

#include "DumpToFile.hpp"

using namespace RRStrategy;

error_code 
DumpToFile::Process () {

	std::stringstream fname;
	const char* uid = Attribute ("UID");

	printf ("Dumping ...\n");

	if (uid == 0  ||  uid == "")
		uid = "unspecified";

	fname << uid << "_config.xml";
	DumpConfig (fname.str().c_str());

	fname.str("");

	fname << uid << "_data.h5";

	while (!m_cplx.empty()) {

		std::stringstream fn;
		
		fn << fname << m_cplx.begin()->first << "_cplx_" << uid << ".h5";
		std::cout << "dumping " << fn.str().c_str() << std::endl;
		m_cplx.begin()->second->Dump(fn.str().c_str());
		m_cplx.erase(m_cplx.begin());

	}

	//m_raw.Dump    (fname.str().c_str());
	//m_helper.Dump (fname.str().c_str());
	//m_pixel.Dump  (fname.str().c_str());

	printf ("... done\n");

	return RRSModule::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new DumpToFile;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

