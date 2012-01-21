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
#include "ReconContext.hpp"
#include "Matrix.hpp"
#include "DataBase.hpp"

using namespace RRServer;

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
ReconContext::SetCXFL          (const std::string name, const cxfl_data& r) {
	DataBase::Instance()->SetCXFL(name, r);
}
		

void
ReconContext::SetCXFL          (const std::string name, Matrix<cxfl>& r) {
	DataBase::Instance()->SetCXFL(name, r);
}
		

void
ReconContext::GetCXFL          (const std::string name, cxfl_data& r) {
	DataBase::Instance()->GetCXFL(name, r);
}
		

void
ReconContext::GetCXFL          (const std::string name, Matrix<cxfl>& r) {
	DataBase::Instance()->GetCXFL(name, r);
}
		

void
ReconContext::SetCXDB          (const std::string name, const cxdb_data& r) {
	DataBase::Instance()->SetCXDB(name, r);
}
		

void
ReconContext::SetCXDB          (const std::string name, Matrix<cxdb>& r) {
	DataBase::Instance()->SetCXDB(name, r);
}
		

void
ReconContext::GetCXDB          (const std::string name, cxdb_data& r) {
	DataBase::Instance()->GetCXDB(name, r);
}
		

void
ReconContext::GetCXDB          (const std::string name, Matrix<cxdb>& r) {
	DataBase::Instance()->GetCXDB(name, r);
}
		

void
ReconContext::SetRLDB (const std::string name, const rldb_data& r) {
	DataBase::Instance()->SetRLDB(name, r);
}
		

void
ReconContext::SetRLDB (const std::string name, Matrix<double>& r) {
	DataBase::Instance()->SetRLDB(name, r);
}
		

void
ReconContext::GetRLDB (const std::string name, Matrix<double>& r) {
	DataBase::Instance()->GetRLDB(name, r);
}
		

void
ReconContext::GetRLDB (const std::string name, rldb_data& r) {
	DataBase::Instance()->GetRLDB(name, r);
}
		

void
ReconContext::SetRLFL (const std::string name, const rlfl_data& r) {
	DataBase::Instance()->SetRLFL(name, r);
}
		

void
ReconContext::SetRLFL (const std::string name, Matrix<float>& r) {
	DataBase::Instance()->SetRLFL(name, r);
}
		

void
ReconContext::GetRLFL (const std::string name, Matrix<float>& r) {
	DataBase::Instance()->GetRLFL(name, r);
}
		

void
ReconContext::GetRLFL (const std::string name, rlfl_data& r) {
	DataBase::Instance()->GetRLFL(name, r);
}
		

void
ReconContext::SetSHRT (const std::string name, const shrt_data& r) {
	DataBase::Instance()->SetSHRT(name, r);
}
		

void
ReconContext::SetSHRT (const std::string name, Matrix<short>& r) {
	DataBase::Instance()->SetSHRT(name, r);
}
		

void
ReconContext::GetSHRT (const std::string name, Matrix<short>& r) {
	DataBase::Instance()->GetSHRT(name, r);
}
		

void
ReconContext::GetSHRT (const std::string name, shrt_data& r) {
	DataBase::Instance()->GetSHRT(name, r);
}
		

void
ReconContext::SetLONG (const std::string name, const long_data& r) {
	DataBase::Instance()->SetLONG(name, r);
}
		

void
ReconContext::SetLONG (const std::string name, Matrix<long>& r) {
	DataBase::Instance()->SetLONG(name, r);
}
		

void
ReconContext::GetLONG (const std::string name, Matrix<long>& r) {
	DataBase::Instance()->GetLONG(name, r);
}
		

void
ReconContext::GetLONG (const std::string name, long_data& r) {
	DataBase::Instance()->GetLONG(name, r);
}
		

void 
ReconContext::Name (const char* name) { 
	m_strategy->Name(name);
}
		
		
const char* 
ReconContext::Name () {return m_strategy->Name();}




