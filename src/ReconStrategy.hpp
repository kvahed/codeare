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

/**
 * @brief Reconstruction / Manipulation modules
 */
namespace RRStrategy {

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
		 * @brief       Default virtul destructor
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
		 * @brief       Mandatory implementation of preparation hooks
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
		Finalise        () { return OK; };
	
	
		/**
		 * @brief       Reference to Database singleton
		 */
		DataBase&
		DB              () const {
			return *(DataBase::Instance());
		}


		/**
		 * @brief       Add a matrix to database map
		 *
		 * @param  name Name in map
		 * @param  p    Smart pointer to the matrix
		 * @return      Reference to matrix
		 */
		template <class T> Matrix<T>& 
		AddMatrix         (const string name, Ptr< Matrix<T> > p) const {
			return DataBase::Instance()->AddMatrix(name, p);
		}
		
		
		/**
		 * @brief       Get reference to complex single matrix by name from database 
		 *              @see DataBase::GetCXFL(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to data 
		 */
		Matrix<cxfl>& 
		GetCXFL         (const string name) const {
			return DataBase::Instance()->GetCXFL(name);
		}
		
		
		/**
		 * @brief       Clear database of complex single matrix by name
		 *              @see DataBase::FreeCXFL(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		bool 
		FreeCXFL        (const string name) const {
			return DataBase::Instance()->FreeCXFL(name);
		}


		/**
		 * @brief       Get reference to complex double matrix by name from database
		 *              @see DataBase::GetCXDB(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		Matrix<cxdb>& 
		GetCXDB         (const string name) const {
			return DataBase::Instance()->GetCXDB(name);
		}
		
		
		/**
		 * @brief       Clear database of complex double matrix by name
		 *              @see DataBase::FreeCXDB(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		bool 
		FreeCXDB        (const string name) const {
			return DataBase::Instance()->FreeCXDB(name);
		}


		/**
		 * @brief       Get reference to single matrix by name from database
		 *              @see DataBase::GetRLFL(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to matrix if existent
		 */
		Matrix<float>& 
		GetRLFL         (const string name) const {
			return DataBase::Instance()->GetRLFL(name);
		}
		
		
		/**
		 * @brief       Clear database of single matrix by name
		 *              @see DataBase::FreeRLFL(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to matrix if existent
		 */
		bool 
		FreeRLFL        (const string name) const {
			return DataBase::Instance()->FreeRLFL(name);
		}


		/**
		 * @brief       Get reference to double matrix by name from database
		 *              @see DataBase::GetRLDB(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to matrix if existent
		 */
		Matrix<double>& 
		GetRLDB         (const string name) const {
			return DataBase::Instance()->GetRLDB(name);
		}
		
		
		/**
		 * @brief       Clear database of double matrix by name
		 *              @see DataBase::FreeRLDB(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to matrix if existent
		 */
		bool 
		FreeRLDB        (const string name) const {
			return DataBase::Instance()->FreeRLDB(name);
		}


		/**
		 * @brief       Get reference to short int matrix by name from database
		 *              @see DataBase::GetSHRT(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to matrix if existent
		 */
		Matrix<short>& 
		GetSHRT         (const string name) const {
			return DataBase::Instance()->GetSHRT(name);
		}
		
		
		/**
		 * @brief       Clear database of short int matrix by name
		 *              @see DataBase::FreeSHRT(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to matrix if existent
		 */
		bool 
		FreeSHRT        (const string name) const {
			return DataBase::Instance()->FreeSHRT(name);
		}


		/**
		 * @brief       Get reference to long int matrix by name from database
		 *              @see DataBase::GetLONG(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to matrix if existent
		 */
		Matrix<long>& 
		GetLONG         (const string name) const {
			return DataBase::Instance()->GetLONG(name);
		}
		
		
		/**
		 * @brief       Clear database of long int matrix by name
		 *              @see DataBase::FreeLONG(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to matrix if existent
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
typedef RRStrategy::ReconStrategy* create_t  ();


/**
 *@brief               Dynamic destructor
 */
typedef void           destroy_t (RRStrategy::ReconStrategy*);


