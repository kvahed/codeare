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

error_code 
DumpToFile::ProcessData () {

	printf ("Dumping ... ");

	if (m_have_raw) {
		printf ("raw data.\n");
		m_raw.dump    ((char*) "raw.bin"   );
	}

	if (m_have_helper) {
		printf ("helper data.\n");
		m_helper.dump ((char*) "helper.bin");
	}

	if (m_have_pixel) {
		printf ("pixel data.\n");
		m_pixel.dump  ((char*) "pixel.bin" );
	}

	if (m_labels.length() > 0) {
		printf ("pixel data.\n");

	    ofstream of;

	    of.open ("labels.txt");
		for (int i = 0; i < m_labels.length(); i++)
			of << m_labels[i] << "\n";
	    of.close();

	}

	return RRSModule::OK;

}
