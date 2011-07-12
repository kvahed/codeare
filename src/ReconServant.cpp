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
ReconServant::ReconServant  ()               {

	m_config = new char;

	m_raw.dims.length(INVALID_DIM);
	m_rhelper.dims.length(INVALID_DIM);
	m_helper.dims.length(INVALID_DIM);
	m_kspace.dims.length(INVALID_DIM);
    m_pixel.dims.length(INVALID_DIM);

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

	e = m_contexts.at(s)->Init();

	std::cout << "... done. Processing ... " << std::endl;

	e = m_contexts.at(s)->Process();
	
	std::cout << "... done. Getting processed data..." << std::endl;
	
	m_contexts.at(s)->GetRaw(&m_raw);
	m_contexts.at(s)->GetRHelper(&m_rhelper);
	m_contexts.at(s)->GetHelper(&m_helper);
	m_contexts.at(s)->GetKSpace(&m_kspace);
	m_contexts.at(s)->GetPixel(&m_pixel);
	m_contexts.at(s)->GetConfig(m_config);
	
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

	return new raw_data (m_raw);

}


/**************************************************************************************************/
void
ReconServant::rhelper          (const raw_data& d)   {

	m_contexts.at(0)->SetRHelper(&d);

}


/**************************************************************************************************/
raw_data*
ReconServant::rhelper          ()                    {
	return new raw_data (m_rhelper);
}


/**************************************************************************************************/
void
ReconServant::helper       (const helper_data& d)   {

	m_contexts.at(0)->SetHelper(&d);

}


/**************************************************************************************************/
helper_data*
ReconServant::helper       ()                    {
	return new helper_data (m_helper);
}


/**************************************************************************************************/
void
ReconServant::kspace       (const helper_data& d)   {

	m_contexts.at(0)->SetKSpace(&d);

}


/**************************************************************************************************/
helper_data*
ReconServant::kspace       ()                    {
	return new helper_data (m_kspace);
}


/**************************************************************************************************/
void
ReconServant::pixel        (const pixel_data& d) {

	m_contexts.at(0)->SetPixel(&d);

}


/**************************************************************************************************/
pixel_data*
ReconServant::pixel        ()                    {
	return new pixel_data (m_pixel);
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

