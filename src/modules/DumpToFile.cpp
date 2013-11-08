/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
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
#include "Algos.hpp"

using namespace RRStrategy;
error_code 
DumpToFile::Process () {

	std::stringstream fname;
	const char* uid = Attribute ("UID");
/*
	map<string,Ptr< Matrix<cxfl> > >&   cfm = Workspace::Instance()->CXFLMap();
	map<string,Ptr< Matrix<double> > >& rdm = Workspace::Instance()->RLDBMap();
	map<string,Ptr< Matrix<short> > >&  sim = Workspace::Instance()->SHRTMap();
*/

	printf ("Dumping ...\n");

	if (uid == 0  || uid[0] == '\0')
		uid = "unspecified";

	fname << uid << "_config.xml";
	DumpConfig (fname.str().c_str());

	fname.str("");

	fname << uid << "_data.mat";
/*
#ifdef HAVE_MAT_H
	MATFile* mf = matOpen (fname.str().c_str(), "w");

	if (mf == NULL) {
		printf ("Error creating file %s\n", fname.str().c_str());
		return FILE_ACCESS_FAILED;
	}

	map<string,Ptr< Matrix<cxfl> > >::iterator cit = cfm.begin();
	for (cit = cfm.begin() ; cit != cfm.end(); cit++) {
		cout << "Dumping " << cit->first.c_str() << endl;
		MXDump(*(cit->second), mf, cit->first.c_str());
	}
	
	map<string,Ptr< Matrix<double> > >::iterator rit = rdm.begin();
	for (rit = rdm.begin(); rit != rdm.end(); rit++) {
		cout << "Dumping " <<  rit->first.c_str() << endl;
		MXDump(*(rit->second), mf, rit->first.c_str());
	}
	
	map<string,Ptr< Matrix<short> > >::iterator pit = sim.begin();
	for (pit = sim.begin(); pit != sim.end(); pit++) {
		cout << "Dumping " <<  pit->first.c_str() << endl;
		MXDump(*(pit->second), mf, pit->first.c_str());
	}

	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n",fname.str().c_str());
		return FILE_ACCESS_FAILED;
	}
#endif
*/
	printf ("... done\n");

	return OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new DumpToFile;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

