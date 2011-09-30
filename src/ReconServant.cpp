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

#include "ReconServant.hpp"
#include "ReconContext.hpp"

using namespace RRServer;

ReconServant::ReconServant  () : 
	m_config (new char) {

}


ReconServant::~ReconServant ()               {

	this->Finalise (-1);
	delete m_config;

}


short
ReconServant::Init (const char* name) {

	error_code e = OK;
	
	m_contexts.push_back(new ReconContext(name));

	std::cout << "Configuring " << m_contexts.at(0)->Name() << " ... " << std::endl;

	m_contexts.at(0)->SetConfig (m_config);

	std::cout << "... done. Initialising algorithm ... " << std::endl;

	if ((e = m_contexts.at(0)->Init()) != OK) {
		Finalise(-1);
		std::cout << "... FAILED!!! Bailing out." << std::endl;
		return -1;
	}

	std::cout << "... done. Handling back control to client ... " << std::endl;

	return e;

}


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


error_code
ReconServant::Process  (const short s)       {

	std::cout << "Processing ... " << std::endl;
	error_code e = m_contexts.at(s)->Process();
	std::cout << "... done. " << std::endl;
	
	return e;
	
}


error_code
ReconServant::Prepare  (const short s)       {
	
	std::cout << "Preparing ... " << std::endl;
	error_code e = m_contexts.at(s)->Prepare();
	std::cout << "... done. " << std::endl;
	
	return e;
	
}


void
ReconServant::set_cplx  (const char* name, const cplx_data& c) {

	m_contexts.at(0)->SetCplx(name, c);

}


void
ReconServant::get_cplx (const char* name, cplx_data& c) {

	c.dims.length(INVALID_DIM);
	c.res.length(INVALID_DIM);
	m_contexts.at(0)->GetCplx(name, c);

}


void
ReconServant::set_real       (const char* name, const real_data& r)   {

	m_contexts.at(0)->SetReal(name, r);

}


void
ReconServant::get_real       (const char* name, real_data& r) {

	r.dims.length(INVALID_DIM);
	r.res.length(INVALID_DIM);
	m_contexts.at(0)->GetReal(name, r);

}


void
ReconServant::set_pixel        (const char* name, const pixel_data& p) {

	m_contexts.at(0)->SetPixel(name, p);

}


void
ReconServant::get_pixel        (const char* name, pixel_data& p) {

	p.dims.length(INVALID_DIM);
	p.res.length(INVALID_DIM);
	m_contexts.at(0)->GetPixel(name, p);

}


void 
ReconServant::config       (const char* d)    {

	std::stringstream tmp;

	tmp << d;
	m_config = new char[tmp.str().length() + 1];
	strcpy (m_config, tmp.str().c_str());

}


char* 
ReconServant::config       ()                    {

	std::stringstream tmp;

	tmp << m_config;

	return CORBA::string_dup(tmp.str().c_str());

}

