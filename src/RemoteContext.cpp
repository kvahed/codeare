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

#include "Loader.hpp"
#include "RemoteContext.hpp"
#include "Matrix.hpp"
#include "DataBase.hpp"

using namespace RRServer;

RemoteContext::~RemoteContext () {
	
	std::cout << "Deleting context " << m_strategy->Name() << " ..." << std::endl;

	destroy_t* destroy = (destroy_t*) GetFunction(m_dlib, (char*)"destroy");
	destroy(m_strategy);
	CloseModule (m_dlib);

	std::cout << "... done." << std::endl << "-----------------------------" << std::endl; 
	
}


RemoteContext::RemoteContext (const char* name) {

	m_dlib = LoadModule ((char*)name);
	create_t* create = (create_t*) GetFunction (m_dlib, (char*)"create");
	m_strategy = create();
	m_strategy->Name (name);
	
}



RemoteContext::RemoteContext     () {}
		
		
RRSModule::error_code
RemoteContext::Process          () {
	return m_strategy->Process();
}


RRSModule::error_code
RemoteContext::Init             () {
	return m_strategy->Init();
}


RRSModule::error_code
RemoteContext::Prepare             () {
	return m_strategy->Prepare();
}


RRSModule::error_code
RemoteContext::Finalise     () {
	return m_strategy->Finalise();
}


void
RemoteContext::SetConfig        (const char* cstr) {
	return m_strategy->SetConfig(cstr);
}


void
RemoteContext::ReadConfig       (const char* fname) {
	m_strategy->ReadConfig(fname);
}


void
RemoteContext::SetCXFL          (const std::string name, const cxfl_data& r) {
	DataBase::Instance()->SetCXFL(name, r);
}


void
RemoteContext::GetCXFL          (const std::string name, cxfl_data& r) {
	DataBase::Instance()->GetCXFL(name, r);
}


void
RemoteContext::SetCXDB          (const std::string name, const cxdb_data& r) {
	DataBase::Instance()->SetCXDB(name, r);
}


void
RemoteContext::GetCXDB          (const std::string name, cxdb_data& r) {
	DataBase::Instance()->GetCXDB(name, r);
}


void
RemoteContext::SetRLDB (const std::string name, const rldb_data& r) {
	DataBase::Instance()->SetRLDB(name, r);
}


void
RemoteContext::GetRLDB (const std::string name, rldb_data& r) {
	DataBase::Instance()->GetRLDB(name, r);
}


void
RemoteContext::SetRLFL (const std::string name, const rlfl_data& r) {
	DataBase::Instance()->SetRLFL(name, r);
}


void
RemoteContext::GetRLFL (const std::string name, rlfl_data& r) {
	DataBase::Instance()->GetRLFL(name, r);
}


void
RemoteContext::SetSHRT (const std::string name, const shrt_data& r) {
	DataBase::Instance()->SetSHRT(name, r);
}


void
RemoteContext::GetSHRT (const std::string name, shrt_data& r) {
	DataBase::Instance()->GetSHRT(name, r);
}


void
RemoteContext::SetLONG (const std::string name, const long_data& r) {
	DataBase::Instance()->SetLONG(name, r);
}


void
RemoteContext::GetLONG (const std::string name, long_data& r) {
	DataBase::Instance()->GetLONG(name, r);
}


void 
RemoteContext::Name (const char* name) { 
	m_strategy->Name(name);
}


const char* 
RemoteContext::Name () {return m_strategy->Name();}




