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

#include "Loader.hpp"
#include "ReconContext.hpp"
#include "Matrix.hpp"
#include "DataBase.hpp"

namespace RRStrategy {


ReconContext::~ReconContext () {
	
	std::cout << "Deleting context " << m_strategy->Name() << " ..." << std::endl;

	destroy_t* destroy = (destroy_t*) GetFunction(m_dlib, (char*)"destroy");
	destroy(m_strategy);
	CloseModule (m_dlib);

	std::cout << "... done." << std::endl << "-----------------------------" << std::endl; 
	
}


ReconContext::ReconContext (const char* name) {

	m_dlib = LoadModule ((char*)name);
	create_t* create = (create_t*) GetFunction (m_dlib, (char*)"create");
	m_strategy = create();
	m_strategy->Name (name);
	
}



ReconContext::ReconContext     () {}
		
		
RRSModule::error_code
ReconContext::Process          () {
	return m_strategy->Process();
}


RRSModule::error_code
ReconContext::Init             () {
	return m_strategy->Init();
}


RRSModule::error_code
ReconContext::Prepare             () {
	return m_strategy->Prepare();
}


RRSModule::error_code
ReconContext::Finalise     () {
	return m_strategy->Finalise();
}


void
ReconContext::SetConfig        (const char* cstr) {
	return m_strategy->SetConfig(cstr);
}


void
ReconContext::ReadConfig       (const char* fname) {
	m_strategy->ReadConfig(fname);
}


void 
ReconContext::Name (const char* name) { 
	m_strategy->Name(name);
}


const char* 
ReconContext::Name () {return m_strategy->Name();}


}

