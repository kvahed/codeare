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
#include "Workspace.hpp"

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
		ReconStrategy   () : m_initialised (false) {}
		

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
		Workspace&
		DB              () const {
			return *(Workspace::Instance());
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
			return Workspace::Instance()->AddMatrix(name, p);
		}
		
		
		/***
		 * @brief       Add a matrix to workspace
		 *
		 * @param  name Name in workspace
		 * @param  vt   Data type of matrix elements
		 * @param  dims Dimensions 
		 * @return      Reference to matrix
		 */
		//template <class T> Matrix<T>& 
		//AddMatrix         (const string& name, const data_type dt, std::vector<size_t> dims) const {
		//		return Workspace::Instance()->AddMatrix(name, dt, dims);
		//}
		
		
		/**
		 * @brief       Get reference to complex single matrix by name from database 
		 *              @see Workspace::Get<T>(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to data 
		 */
		template <class T> Matrix<T>& 
		Get            (const string name) const {
			return Workspace::Instance()->Get<T>(name);
		}
		
		
		/**
		 * @brief       Clear database of complex single matrix by name
		 *              @see Workspace::FreeCXFL(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		inline bool 
		Free            (const string name) const {
			return Workspace::Instance()->Free (name);
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


