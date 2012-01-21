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

#ifndef __RECON_STRATEGY_HPP__
#define __RECON_STRATEGY_HPP__

#include "Matrix.hpp"
#include "Configurable.hpp"
#include "DataBase.hpp"

#include "DllExport.h"

#ifdef __WIN32__ 
    #include "RRSModule.h"
#else
    #include "RRSModule.hh"
#endif

#include <cstdlib>
#include <complex>
#include <stdint.h>

using namespace RRSModule;
using namespace std;

namespace RRServer {

	/**
	 * @brief Strategy for reconstruction strategies
	 *        Derive hereof to expand the reconstruction toolbox
	 *
	 */
	class ReconStrategy : public Configurable {
		
		
	public:
		

		/**
		 * @brief       Default constructor
		 */ 
		ReconStrategy   () {}
		

		/**
		 * @brief       Default destructor
		 */ 
		virtual
		~ReconStrategy  () {}
		
		
		/**
		 * @brief       Mandatory implementation of actual data procession
		 *
		 * @return      Success
		 */ 
		virtual error_code
		Process         () = 0;
		

		/**
		 * @brief       Mandatory implementation of initialiser
		 *
		 * @return      Success
		 */ 
		virtual error_code
		Init            () = 0;
		

		/**
		 * @brief       Additional preparation?
		 *
		 * @return      Success
		 */ 
		virtual error_code
		Prepare         () { 
			return RRSModule::OK; 
		}
		

		/**
		 * @brief       Attach a name to the algorithm
		 *
		 * @param  name Name
		 */
		void 
		Name            (const char* name) { 
			m_name = string(name);
		}


		/**
		 * @brief       Get given name
		 *
		 * @return      Name
		 */
		const char* 
		Name            () {
			return m_name.c_str();
		}


		/**
		 * @brief       Clean up
		 *
		 * @return      Success
		 */ 
		virtual error_code 
		Finalise        () {};
	
	
		DataBase&
		DB              () const {
			return *(DataBase::Instance());
		}


		/**
		 * @brief       Add a complex matrix to complex single matrix database
		 *
		 * @param  name Name
		 * @param  m    The added matrix
		 * @return      Success
		 */
		Matrix<cplx>& 
		AddCXFL         (const string name, Ptr< Matrix<cplx> > m) const {
			return DataBase::Instance()->AddCXFL(name, m);
		}
		
		
		/**
		 * @brief       Get reference to complex single data by name from database
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		Matrix<cplx>& 
		GetCXFL         (const string name) const {
			return DataBase::Instance()->GetCXFL(name);
		}
		
		
		/**
		 * @brief       Clear database of complex single data by name
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		bool 
		FreeCXFL        (const string name) const {
			return DataBase::Instance()->FreeCXFL(name);
		}


		/**
		 * @brief       Add a complex matrix to complex double matrix database
		 *
		 * @param  name Name
		 * @param  m    The added matrix
		 * @return      Success
		 */
		Matrix<cxdb>& 
		AddCXDB         (const string name, Ptr< Matrix<cxdb> > m) const {
			return DataBase::Instance()->AddCXDB(name, m);
		}
		
		
		/**
		 * @brief       Get reference to complex double data by name from database
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		Matrix<cxdb>& 
		GetCXDB         (const string name) const {
			return DataBase::Instance()->GetCXDB(name);
		}
		
		
		/**
		 * @brief       Clear database of complex double data by name
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		bool 
		FreeCXDB        (const string name) const {
			return DataBase::Instance()->FreeCXDB(name);
		}


		/**
		 * @brief       Add a complex matrix to complex single matrix database
		 *
		 * @param  name Name
		 * @param  m    The added matrix
		 * @return      Success
		 */
		Matrix<float>& 
		AddRLFL         (const string name, Ptr< Matrix<float> > m) const {
			return DataBase::Instance()->AddRLFL(name, m);
		}
		
		
		/**
		 * @brief       Get reference to complex single data by name from database
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		Matrix<float>& 
		GetRLFL         (const string name) const {
			return DataBase::Instance()->GetRLFL(name);
		}
		
		
		/**
		 * @brief       Clear database of complex single data by name
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		bool 
		FreeRLFL        (const string name) const {
			return DataBase::Instance()->FreeRLFL(name);
		}


		/**
		 * @brief       Add a complex matrix to complex double matrix database
		 *
		 * @param  name Name
		 * @param  m    The added matrix
		 * @return      Success
		 */
		Matrix<double>& 
		AddRLDB         (const string name, Ptr< Matrix<double> > m) const {
			return DataBase::Instance()->AddRLDB(name, m);
		}
		
		
		/**
		 * @brief       Get reference to complex double data by name from database
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		Matrix<double>& 
		GetRLDB         (const string name) const {
			return DataBase::Instance()->GetRLDB(name);
		}
		
		
		/**
		 * @brief       Clear database of complex double data by name
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		bool 
		FreeRLDB        (const string name) const {
			return DataBase::Instance()->FreeRLDB(name);
		}


		/**
		 * @brief       Add a complex matrix to complex double matrix database
		 *
		 * @param  name Name
		 * @param  m    The added matrix
		 * @return      Success
		 */
		Matrix<short>& 
		AddSHRT         (const string name, Ptr< Matrix<short> > m) const {
			return DataBase::Instance()->AddSHRT(name, m);
		}
		
		
		/**
		 * @brief       Get reference to complex double data by name from database
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		Matrix<short>& 
		GetSHRT         (const string name) const {
			return DataBase::Instance()->GetSHRT(name);
		}
		
		
		/**
		 * @brief       Clear database of complex double data by name
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		bool 
		FreeSHRT        (const string name) const {
			return DataBase::Instance()->FreeSHRT(name);
		}


		/**
		 * @brief       Add a complex matrix to complex double matrix database
		 *
		 * @param  name Name
		 * @param  m    The added matrix
		 * @return      Success
		 */
		Matrix<long>& 
		AddLONG         (const string name, Ptr< Matrix<long> > m) const {
			return DataBase::Instance()->AddLONG(name, m);
		}
		
		
		/**
		 * @brief       Get reference to complex double data by name from database
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		Matrix<long>& 
		GetLONG         (const string name) const {
			return DataBase::Instance()->GetLONG(name);
		}
		
		
		/**
		 * @brief       Clear database of complex double data by name
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		bool 
		FreeLONG        (const string name) const {
			return DataBase::Instance()->FreeLONG(name);
		}


	protected:
		
		string    m_name;         /*!< @brief Name                        */
		bool      m_initialised;  /*!< @brief Reco is initialised         */
		
	};
	
}
#endif /* __RECON_STRATEGY_HPP__ */
	


/**
 * @brief              Dynamic constructor
 */
typedef RRServer::ReconStrategy* create_t  ();


/**
 *@brief               Dynamic destructor
 */
typedef void           destroy_t (RRServer::ReconStrategy*);


