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
#include "SimpleTimer.hpp"

#include "DllExport.h"


#include <cstdlib>
#include <complex>

class Workspace;

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
		ReconStrategy   () : m_initialised (false), global (0) {}
		

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
		virtual codeare::error_code
		Process         () = 0;
		

		/**
		 * @brief       Mandatory implementation of initialiser
		 *
		 * @return      Success
		 */ 
		virtual codeare::error_code
		Init            () = 0;
		

		/**
		 * @brief       Mandatory implementation of preparation hooks
		 *
		 * @return      Success
		 */ 
		virtual codeare::error_code
		Prepare         () { 
			return codeare::OK; 
		}
		

		/**
		 * @brief       Attach a name to the algorithm
		 *
		 * @param  name Name
		 */
		void 
		Name            (const char* name) { 
			m_name = std::string(name);
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
		virtual codeare::error_code 
		Finalise        () { return codeare::OK; };
	
	
		/**
		 * @brief       Reference to Database singleton
		 */
		Workspace&
		DB              () const {
			return *global;
		}


		/**
		 * @brief       Add a matrix to database map
		 *
		 * @param  name Name in map
		 * @param  p    Smart pointer to the matrix
		 * @return      Reference to matrix
		 */
		template <class T> Matrix<T>& 
		AddMatrix         (const std::string& name, boost::shared_ptr< Matrix<T> > p) const {
			return global->AddMatrix(name, p);
		}


		/**
		 * @brief       Add a matrix to database map
		 *
		 * @param  name Name in map
		 * @param  p    Smart pointer to the matrix
		 * @return      Reference to matrix
		 */
		template <class T> Matrix<T>& 
		AddMatrix         (const std::string& name) const {
			return global->AddMatrix<T>(name);
		}


		/**
		 * @brief       Add a matrix to database map
		 *
		 * @param  name Name in map
		 * @param  p    Smart pointer to the matrix
		 * @return      Reference to matrix
		 */
		template <class T> Matrix<T>& 
		AddMatrix         (const char* name) const {
			return global->AddMatrix<T>(std::string(name));
		}


		/**
		 * @brief       Get reference to complex single matrix by name from database
		 *              @see Workspace::Get<T>(const string)
		 *
		 * @param  name Name
		 * @return      Reference to data
		 */
		template <class T> Matrix<T>&
		Get            (const std::string& name) const {
			return global->Get<T>(name);
		}
		
		
		/**
		 * @brief       Clear database of complex single matrix by name
		 *              @see Workspace::FreeCXFL(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		inline bool 
		Free            (const std::string& name) const {
			return global->Free (name);
		}


		/**
		 * @brief       Clear database of complex single matrix by name
		 *              @see Workspace::FreeCXFL(const string)
		 * 
		 * @param  name Name
		 * @return      Reference to data if existent
		 */
		inline void 
		WSpace         (Workspace* ws) {
			global = ws;
		}


	protected:
		
		std::string    m_name;         /*!< @brief Name                        */
		bool           m_initialised;  /*!< @brief Reco is initialised         */
		ReconStrategy* _successor;

        Workspace*     global;
		
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


