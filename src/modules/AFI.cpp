/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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

#include "AFI.hpp"

using namespace RRStrategy;



AFI::AFI() : m_use_real(true), 
			 m_retain_phase (true) {



}


AFI::~AFI() {

	this->Finalise();

}


RRSModule::error_code
AFI::Init () {

	Attribute ("use_real",     &m_use_real);
	Attribute ("retain_phase", &m_retain_phase);

	m_initialised = true;
	return RRSModule::OK;

}


RRSModule::error_code
AFI::Finalise () {

	return ReconStrategy::Finalise();
	
}


RRSModule::error_code
AFI::Process     () { 

	Matrix<cxfl>& afid = GetCXFL("meas");
	
	Matrix<float> ph1 (afid.Dim(0), afid.Dim(1), afid.Dim(2));
	PhasePreset (afid, m_use_real, ph1);
	
	return RRSModule::OK;
	
}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {

    return new AFI;

}


extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {

	delete p;

}
