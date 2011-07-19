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

#include "ReconServant.hpp"

using namespace RRServer;

/**************************************************************************************************/
ReconServant::ReconServant  () : 
	m_config (new char) {

}


/**************************************************************************************************/
ReconServant::~ReconServant ()               {

	this->Finalise (-1);
	delete m_config;

}


/**************************************************************************************************/
short
ReconServant::Init (const char* name) {

	error_code e = OK;
	
	m_contexts.push_back(new ReconContext(name));

	return e;

}


/**************************************************************************************************/
error_code
ReconServant::Finalise (const short s) {

	error_code e = OK;

	if (s == -1)
		for (int i = 0; i < m_contexts.size(); i++) {
			delete m_contexts.back();
			m_contexts.pop_back();
		}
	else {
		
		delete m_contexts.at(s);
		m_contexts.erase (m_contexts.begin()+s);
		
	}

	return e;

}


/**************************************************************************************************/
error_code
ReconServant::Process  (const short s)       {

	error_code e = OK;
	
	std::cout << "Configuring " << m_contexts.at(s)->Name() << " ... " << std::endl;

	m_contexts.at(s)->SetConfig(m_config);

	std::cout << "... done. Initialising algorithm ... " << std::endl;

	if ((e = m_contexts.at(s)->Init()) != OK) {
		std::cout << "... FAILED!!! Bailing out." << std::endl;
		return e;
	}

	std::cout << "... done. Processing ... " << std::endl;

	e = m_contexts.at(s)->Process();
	
	std::cout << "... done. " << std::endl;
	
	return e;
	
}


/**************************************************************************************************/
void
ReconServant::raw          (const raw_data& d)   {
	
	m_contexts.at(0)->SetRaw(&d);

}


/**************************************************************************************************/
raw_data*
ReconServant::raw          ()                    {

	raw_data tmp;      
	tmp.dims.length(INVALID_DIM);
	m_contexts.at(0)->GetRaw(&tmp);
	return new raw_data (tmp);

}


/**************************************************************************************************/
void
ReconServant::rhelper          (const raw_data& d)   {

	m_contexts.at(0)->SetRHelper(&d);

}


/**************************************************************************************************/
raw_data*
ReconServant::rhelper          ()                    {

	raw_data tmp;      
	tmp.dims.length(INVALID_DIM);
	m_contexts.at(0)->GetRHelper(&tmp);
	return new raw_data (tmp);

}


/**************************************************************************************************/
void
ReconServant::helper       (const helper_data& d)   {

	m_contexts.at(0)->SetHelper(&d);

}


/**************************************************************************************************/
helper_data*
ReconServant::helper       ()                    {

	helper_data tmp;      
	tmp.dims.length(INVALID_DIM);
	m_contexts.at(0)->GetHelper(&tmp);
	return new helper_data (tmp);

}


/**************************************************************************************************/
void
ReconServant::kspace       (const helper_data& d)   {

	m_contexts.at(0)->SetKSpace(&d);

}


/**************************************************************************************************/
helper_data*
ReconServant::kspace       ()                    {

	helper_data tmp;      
	tmp.dims.length(INVALID_DIM);
	m_contexts.at(0)->GetKSpace(&tmp);
	return new helper_data (tmp);

}


/**************************************************************************************************/
void
ReconServant::pixel        (const pixel_data& d) {

	m_contexts.at(0)->SetPixel(&d);

}


/**************************************************************************************************/
pixel_data*
ReconServant::pixel        ()                    {

	pixel_data tmp;      
	tmp.dims.length(INVALID_DIM);
	m_contexts.at(0)->GetPixel(&tmp);
	return new pixel_data (tmp);

}


/**************************************************************************************************/
void 
ReconServant::config       (const char* d)    {

	std::stringstream tmp;

	tmp << d;
	m_config = new char[tmp.str().length()];
	strcpy (m_config, tmp.str().c_str());

}


/**************************************************************************************************/
char* 
ReconServant::config       ()                    {

	std::stringstream tmp;

	tmp << m_config;

	return CORBA::string_dup(tmp.str().c_str());

}

