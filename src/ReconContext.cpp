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

#include "Loader.hpp"
#include "ReconContext.hpp"
#include "Matrix.hpp"

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
ReconContext::Finalise     () {
	return m_strategy->Finalise();
}
		

void
ReconContext::SetConfig        (const char* cstr) {
	m_strategy->SetConfig(cstr);
}
		
		
void
ReconContext::ReadConfig       (const char* fname) {
	m_strategy->ReadConfig(fname);
}
		

void
ReconContext::SetRaw          (const raw_data* r) {
	m_strategy->SetRaw(r);
}
		

void
ReconContext::SetRaw          (const Matrix<raw>* r) {
	m_strategy->SetRaw(r);
}
		

void
ReconContext::GetRaw          (raw_data* r) {
	m_strategy->GetRaw(r);
}
		

void
ReconContext::GetRaw          (Matrix<raw>* r) {
	m_strategy->GetRaw(r);
}
		

void
ReconContext::SetRHelper (const raw_data* r) {
	m_strategy->SetRHelper(r);
}
		

void
ReconContext::SetRHelper (const Matrix<raw>* r) {
	m_strategy->SetRHelper(r);
}
		

void
ReconContext::GetRHelper (raw_data* r) {
	m_strategy->GetRHelper(r);
}
		

void
ReconContext::GetRHelper (Matrix<raw>* r) {
	m_strategy->GetRHelper(r);
}
		

void
ReconContext::SetHelper (const helper_data* r) {
	m_strategy->SetHelper(r);
}
		

void
ReconContext::SetHelper (const Matrix<double>* r) {
	m_strategy->SetHelper(r);
}
		

void
ReconContext::GetHelper (Matrix<double>* r) {
	m_strategy->GetHelper(r);
}
		

void
ReconContext::GetHelper (helper_data* r) {
	m_strategy->GetHelper(r);
}
		

void
ReconContext::SetKSpace (const helper_data* r) {
	m_strategy->SetKSpace(r);
}
		

void
ReconContext::SetKSpace (const Matrix<double>* r) {
	m_strategy->SetKSpace(r);
}
		

void
ReconContext::GetKSpace (Matrix<double>* r) {
	m_strategy->GetKSpace(r);
}
		

void
ReconContext::GetKSpace (helper_data* r) {
	m_strategy->GetKSpace(r);
}
		

void
ReconContext::SetPixel (const pixel_data* r) {
	m_strategy->SetPixel(r);
}
		

void
ReconContext::SetPixel (const Matrix<short>* r) {
	m_strategy->SetPixel(r);
}
		

void
ReconContext::GetPixel (Matrix<short>* r) {
	m_strategy->GetPixel(r);
}
		

void
ReconContext::GetPixel (pixel_data* r) {
	m_strategy->GetPixel(r);
}
		

void 
ReconContext::Name (const char* name) { 
	m_strategy->Name(name);
}
		
		
const char* 
ReconContext::Name () {return m_strategy->Name();}




