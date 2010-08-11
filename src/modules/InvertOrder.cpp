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

#include "InvertOrder.h"

RRSModule::error_code
InvertOrder::Process () {

	/*	RRSModule::floats tmpabs (m_raw.dabs);
	RRSModule::floats tmparg (m_raw.darg);
	
	int     len = tmpabs.length();
	
	for (int i = 0; i < len; i++) {
		tmpabs[i] = m_raw.dabs [len-1-i];
		tmparg[i] = m_raw.darg [len-1-i];
	}
	
	m_raw.dabs = tmpabs;
	m_raw.darg = tmparg;*/
	
	return RRSModule::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new InvertOrder;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}
