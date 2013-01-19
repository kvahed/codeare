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
#include "Workspace.hpp"

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
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];

		for (int i = 0; i < INVALID_DIM; i++) {
			mdims[i] = c.dims[i];
			mress[i] = c.res[i];
		}

		Matrix<cxfl> pm (mdims, mress);

		memcpy (&pm[0], &c.vals[0], pm.Size() * sizeof(cxfl));

		Workspace::Instance()->SetMatrix(name, pm);
		
	}
	
	
	void
	ReconServant::get_cxfl (const char* name, cxfl_data& c) {
		
		Matrix<cxfl> tmp = Workspace::Instance()->Get<cxfl> (name);

		c.dims.length(INVALID_DIM);
		c.res.length (INVALID_DIM);

		for (int j = 0; j < INVALID_DIM; j++) {
			c.dims[j] = tmp.Dim(j);
			c.res[j]  = tmp.Res(j);
		}

		c.vals.length(2 * tmp.Size());

		memcpy (&c.vals[0], &tmp[0], tmp.Size() * sizeof(cxfl));
		
	}
	

	void
	ReconServant::set_cxdb  (const char* name, const cxdb_data& c) {
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];

		for (int i = 0; i < INVALID_DIM; i++) {
			mdims[i] = c.dims[i];
			mress[i] = c.res[i];
		}

		Matrix<cxdb> pm (mdims, mress);

		memcpy (&pm[0], &c.vals[0], pm.Size() * sizeof(cxdb));

		Workspace::Instance()->SetMatrix(name, pm);

	}


	void
	ReconServant::get_cxdb (const char* name, cxdb_data& c) {

		Matrix<cxdb> tmp = Workspace::Instance()->Get<cxdb> (name);

		c.dims.length(INVALID_DIM);
		c.res.length (INVALID_DIM);

		for (int j = 0; j < INVALID_DIM; j++) {
			c.dims[j] = tmp.Dim(j);
			c.res[j]  = tmp.Res(j);
		}

		c.vals.length(2 * tmp.Size());

		memcpy (&c.vals[0], &tmp[0], tmp.Size() * sizeof(cxdb));

	}


	void
	ReconServant::set_rlfl  (const char* name, const rlfl_data& c) {
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];

		for (int i = 0; i < INVALID_DIM; i++) {
			mdims[i] = c.dims[i];
			mress[i] = c.res[i];
		}

		Matrix<float> pm (mdims, mress);

		memcpy (&pm[0], &c.vals[0], pm.Size() * sizeof(float));

		Workspace::Instance()->SetMatrix(name, pm);
		
	}
	
	
	void
	ReconServant::get_rlfl (const char* name, rlfl_data& c) {
		
		Matrix<float> tmp = Workspace::Instance()->Get<float> (name);

		c.dims.length(INVALID_DIM);
		c.res.length (INVALID_DIM);

		for (int j = 0; j < INVALID_DIM; j++) {
			c.dims[j] = tmp.Dim(j);
			c.res[j]  = tmp.Res(j);
		}

		c.vals.length(tmp.Size());

		memcpy (&c.vals[0], &tmp[0], tmp.Size() * sizeof(float));
		
	}
	

	void
	ReconServant::set_rldb  (const char* name, const rldb_data& c) {
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];

		for (int i = 0; i < INVALID_DIM; i++) {
			mdims[i] = c.dims[i];
			mress[i] = c.res[i];
		}

		Matrix<double> pm (mdims, mress);

		memcpy (&pm[0], &c.vals[0], pm.Size() * sizeof(double));

		Workspace::Instance()->SetMatrix(name, pm);

	}


	void
	ReconServant::get_rldb (const char* name, rldb_data& c) {

		Matrix<double> tmp = Workspace::Instance()->Get<double> (name);

		c.dims.length(INVALID_DIM);
		c.res.length (INVALID_DIM);

		for (int j = 0; j < INVALID_DIM; j++) {
			c.dims[j] = tmp.Dim(j);
			c.res[j]  = tmp.Res(j);
		}

		c.vals.length(tmp.Size());

		memcpy (&c.vals[0], &tmp[0], tmp.Size() * sizeof(cxdb));

	}


	void
	ReconServant::set_shrt  (const char* name, const shrt_data& c) {
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];

		for (int i = 0; i < INVALID_DIM; i++) {
			mdims[i] = c.dims[i];
			mress[i] = c.res[i];
		}

		Matrix<short> pm (mdims, mress);

		memcpy (&pm[0], &c.vals[0], pm.Size() * sizeof(short));

		Workspace::Instance()->SetMatrix(name, pm);

	}


	void
	ReconServant::get_shrt (const char* name, shrt_data& c) {

		Matrix<short> tmp = Workspace::Instance()->Get<short> (name);

		c.dims.length(INVALID_DIM);
		c.res.length (INVALID_DIM);

		for (int j = 0; j < INVALID_DIM; j++) {
			c.dims[j] = tmp.Dim(j);
			c.res[j]  = tmp.Res(j);
		}

		c.vals.length(tmp.Size());

		memcpy (&c.vals[0], &tmp[0], tmp.Size() * sizeof(cxdb));

	}
	

	void
	ReconServant::set_long  (const char* name, const long_data& c) {
		
		size_t mdims [INVALID_DIM];
		float  mress [INVALID_DIM];

		for (int i = 0; i < INVALID_DIM; i++) {
			mdims[i] = c.dims[i];
			mress[i] = c.res[i];
		}

		Matrix<long> pm (mdims, mress);

		memcpy (&pm[0], &c.vals[0], pm.Size() * sizeof(long));

		Workspace::Instance()->SetMatrix(name, pm);

	}


	void
	ReconServant::get_long (const char* name, long_data& c) {

		Matrix<long> tmp = Workspace::Instance()->Get<long> (name);

		c.dims.length(INVALID_DIM);
		c.res.length (INVALID_DIM);

		for (int j = 0; j < INVALID_DIM; j++) {
			c.dims[j] = tmp.Dim(j);
			c.res[j]  = tmp.Res(j);
		}

		c.vals.length(2 * tmp.Size());

		memcpy (&c.vals[0], &tmp[0], tmp.Size() * sizeof(long));

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
