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

#include "ReconServant.hpp"

using namespace RRStrategy;

namespace RRServer {

	ReconServant::ReconServant  () {}
	
	
	ReconServant::~ReconServant () {}
	
	
	error_code
	ReconServant::CleanUp () {
		
		return FunctorContainer::CleanUp();
		
	}
	
	
	error_code
	ReconServant::Init (const char* name) {
		
		return FunctorContainer::Init (name);
		
	}
	
	
	error_code
	ReconServant::Finalise (const char* name) {
		
		return FunctorContainer::Finalise (name);
		
	}
	
	
	error_code
	ReconServant::Process  (const char* name)       {
		
		return FunctorContainer::Process (name);
		
	}
	
	
	error_code
	ReconServant::Prepare  (const char* name)       {
	
		return FunctorContainer::Prepare (name);
		
	}
	
	
	void
	ReconServant::set_cxfl  (const char* name, const cxfl_data& c) {
		
		Workspace::Instance()->SetMatrix(name, c);
		
	}
	
	
	void
	ReconServant::get_cxfl (const char* name, cxfl_data& c) {
		
		c.dims.length(INVALID_DIM);
		c.res.length (INVALID_DIM);
		Workspace::Instance()->GetMatrix(name, c);
		
	}
	

	void
	ReconServant::set_cxdb  (const char* name, const cxdb_data& c) {
		
		Workspace::Instance()->SetMatrix(name, c);
		
	}
	
	
	void
	ReconServant::get_cxdb (const char* name, cxdb_data& c) {
		
		c.dims.length(INVALID_DIM);
		c.res.length (INVALID_DIM);
		Workspace::Instance()->GetMatrix(name, c);
		
	}
	
	
	void
	ReconServant::set_rldb       (const char* name, const rldb_data& r)   {
		
		Workspace::Instance()->SetMatrix(name, r);
		
	}
	
	
	void
	ReconServant::get_rldb       (const char* name, rldb_data& r) {
		
		r.dims.length(INVALID_DIM);
		r.res.length (INVALID_DIM);
		Workspace::Instance()->GetMatrix(name, r);
		
	}
	
	
	void
	ReconServant::set_rlfl       (const char* name, const rlfl_data& r)   {
		
		Workspace::Instance()->SetMatrix(name, r);
		
	}
	
	
	void
	ReconServant::get_rlfl       (const char* name, rlfl_data& r) {
		
		r.dims.length(INVALID_DIM);
		r.res.length(INVALID_DIM);
		Workspace::Instance()->GetMatrix(name, r);
		
	}
	
	
	void
	ReconServant::set_shrt        (const char* name, const shrt_data& p) {
		
		Workspace::Instance()->SetMatrix(name, p);
		
	}
	
	
	void
	ReconServant::get_shrt        (const char* name, shrt_data& p) {
		
		p.dims.length(INVALID_DIM);
		p.res.length(INVALID_DIM);
		Workspace::Instance()->GetMatrix(name, p);
		
	}
	
	
	void
	ReconServant::set_long        (const char* name, const long_data& p) {
		
		Workspace::Instance()->SetMatrix(name, p);
		
	}
	
	
	void
	ReconServant::get_long        (const char* name, long_data& p) {
		
		p.dims.length(INVALID_DIM);
		p.res.length(INVALID_DIM);
		Workspace::Instance()->GetMatrix(name, p);
		
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
	
};
