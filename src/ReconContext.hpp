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

#ifndef __RECONCONTEXT_HPP__
#define __RECONCONTEXT_HPP__

#include "ReconStrategy.hpp"
#include "Connector.hpp"
#include <vector>

namespace RRServer {


	/**
	 * @brief Context of a reconstruction method. Abstraction layer to algorithm backends. 
	 */
	class ReconContext : public Connector {
		
		
		
	public:
		
		
		/**
		 * @brief        Unload library and destruct.
		 */ 
		~ReconContext    ();
		
		
		/**
		 * @brief        Construct with a strategy name.
		 *               Loads and initialises algorithm. 
		 *               Needs a present config [@see ReconServant::config(const char*)]).
		 *
		 * @name         Name of algorithm
		 */
		ReconContext     (const char* name);
		
		
		/**
		 * @brief        Direct access pointer to underlying algorithm.
		 *
		 * @return       Algorithm
		 */
		inline ReconStrategy*
		Strategy         ();
		

		/**
		 * @brief        Process. @see ReconStrategy::Process()
		 *
		 * @return       Success
		 */
		RRSModule::error_code
		Process          ();
		
		
		/**
		 * @brief        Initialise. @see ReconStrategy::Init()
		 *
		 * @return       Success
		 */
		RRSModule::error_code
		Init             ();
		
		
		/**
		 * @brief        Prepare. @see ReconStrategy::Prepare()
		 *
		 * @return       Success
		 */
		RRSModule::error_code
		Prepare          ();
		
		
		/**
		 * @brief        Finalise. @see ReconStrategy::Finalise()
		 *
		 * @return       Success
		 */
		RRSModule::error_code
		Finalise     ();
		

		/**
		 * @brief        @see Reconstrategy::SetConfig(const char*)
		 *
		 * @param  cstr  Serialised XML configuration
		 */
		void
		SetConfig        (const char* cstr);
		
		
		/**
		 * @brief        @see ReconStrategy::ReadConfig(const char*)
		 *
		 * @param  fname File name
		 */
		void
		ReadConfig       (const char* fname);
		

		/**
		 * @brief        @see ReconStrategy::SetCXFL(const std::string name, const cxfl_data&)
		 *
		 * @param  name  Name
		 * @param  r     Complex data sequence
		 */
		void
		SetCXFL          (const std::string name, const cxfl_data& r);
		

		/**
		 * @brief       @see ReconStrategy::SetCXFL(const std::string name, const Matrix<cxfl>&)
		 *
		 * @param  name  Name
		 * @param  r    Complex data matrix
		 */
		void
		SetCXFL          (const std::string name, Matrix<cxfl>& r);
		

		/**
		 * @brief       @see ReconStrategy::GetCXFL(const std::string name, cxfl_data&)
		 *
		 * @param  name  Name
		 * @param  r    Complex data sequence
		 */
		void
		GetCXFL          (const std::string name, cxfl_data& r);
		

		/**
		 * @brief       @see ReconStrategy::GetCXFL(const std::string name, Matrix<cxfl>*)
		 *
		 * @param  name  Name
		 * @param  r    Complex data matrix
		 */
		void
		GetCXFL          (const std::string name, Matrix<cxfl>& r);
		

		/**
		 * @brief        @see ReconStrategy::SetCXDB(const std::string name, const cxdb_data&)
		 *
		 * @param  name  Name
		 * @param  r     Complex data sequence
		 */
		void
		SetCXDB          (const std::string name, const cxdb_data& r);
		

		/**
		 * @brief       @see ReconStrategy::SetCXDB(const std::string name, const Matrix<cxdb>&)
		 *
		 * @param  name  Name
		 * @param  r    Complex data matrix
		 */
		void
		SetCXDB          (const std::string name, Matrix<cxdb>& r);
		

		/**
		 * @brief       @see ReconStrategy::GetCXDB(const std::string name, cxdb_data&)
		 *
		 * @param  name  Name
		 * @param  r    Complex data sequence
		 */
		void
		GetCXDB          (const std::string name, cxdb_data& r);
		

		/**
		 * @brief       @see ReconStrategy::GetCXDB(const std::string name, Matrix<cxdb>*)
		 *
		 * @param  name  Name
		 * @param  r    Complex data matrix
		 */
		void
		GetCXDB          (const std::string name, Matrix<cxdb>& r);
		

		/**
		 * @brief       @see ReconStrategy::SetRLDB(const std::string name, const rldb_data*)
		 *
		 * @param  name  Name
		 * @param  r    Real data sequence
		 */
		void
		SetRLDB (const std::string name, const rldb_data& r);
		

		/**
		 * @brief       @see ReconStrategy::SetRLDB(const std::string name, const Matrix<double>*)
		 *
		 * @param  name  Name
		 * @param  r    Real data matrix
		 */
		void
		SetRLDB (const std::string name, Matrix<double>& r);
		

		/**
		 * @brief       @see ReconStrategy::SetRLDB(const std::string name, const Matrix<double>*)
		 *
		 * @param  name  Name
		 * @param  r    Real data matrix
		 */
		void
		GetRLDB (const std::string name, Matrix<double>& r);
		

		/**
		 * @brief       @see ReconStrategy::GetRLDB(const std::string name, rldb_data*)
		 *
		 * @param  name  Name
		 * @param  r    Real data sequence
		 */
		void
		GetRLDB (const std::string name, rldb_data& r);
		

		/**
		 * @brief       @see ReconStrategy::SetRLFL(const std::string name, const rlfl_data*)
		 *
		 * @param  name  Name
		 * @param  f    Float data sequence
		 */
		void
		SetRLFL (const std::string name, const rlfl_data& f);
		

		/**
		 * @brief       @see ReconStrategy::SetRLFL(const std::string name, const Matrix<float>*)
		 *
		 * @param  name  Name
		 * @param  f     Float data matrix
		 */
		void
		SetRLFL (const std::string name, Matrix<float>& f);
		

		/**
		 * @brief       @see ReconStrategy::SetRLFL(const std::string name, const Matrix<float>*)
		 *
		 * @param  name  Name
		 * @param  f     Float data matrix
		 */
		void
		GetRLFL (const std::string name, Matrix<float>& f);
		

		/**
		 * @brief       @see ReconStrategy::GetRLFL(const std::string name, rlfl_data*)
		 *
		 * @param  name  Name
		 * @param  f    Float data sequence
		 */
		void
		GetRLFL (const std::string name, rlfl_data& f);
		

		/**
		 * @brief       @see ReconStrategy::SetSHRT(const std::string name, const shrt_data*)
		 *
		 * @param  name  Name
		 * @param  r    Short int data sequence
		 */
		void
		SetSHRT (const std::string name, const shrt_data& r);
		

		/**
		 * @brief       @see ReconStrategy::SetSHRT(const std::string name, const Matrix<short>*)
		 *
		 * @param  name  Name
		 * @param  r    Short int data matrix
		 */
		void
		SetSHRT (const std::string name, Matrix<short>& r);
		

		/**
		 * @brief       @see ReconStrategy::GetSHRT(const std::string name, Matrix<short>*)
		 *
		 * @param  name  Name
		 * @param  r    Short int data matrix
		 */
		void
		GetSHRT (const std::string name, Matrix<short>& r);
		

		/**
		 * @brief       @see ReconStrategy::GetSHRT(const std::string name, shrt_data*)
		 *
		 * @param  name  Name
		 * @param  r    Short int data sequence
		 */
		void
		GetSHRT (const std::string name, shrt_data& r);
		

		/**
		 * @brief       @see ReconStrategy::SetLONG(const std::string name, const long_data*)
		 *
		 * @param  name  Name
		 * @param  r    Long int data sequence
		 */
		void
		SetLONG (const std::string name, const long_data& r);
		

		/**
		 * @brief       @see ReconStrategy::SetLONG(const std::string name, const Matrix<long>*)
		 *
		 * @param  name  Name
		 * @param  r    Long int data matrix
		 */
		void
		SetLONG (const std::string name, Matrix<long>& r);
		

		/**
		 * @brief       @see ReconStrategy::GetLONG(const std::string name, Matrix<long>*)
		 *
		 * @param  name  Name
		 * @param  r    Long int data matrix
		 */
		void
		GetLONG (const std::string name, Matrix<long>& r);
		

		/**
		 * @brief       @see ReconStrategy::GetLONG(const std::string name, long_data*)
		 *
		 * @param  name  Name
		 * @param  r    Long int data sequence
		 */
		void
		GetLONG (const std::string name, long_data& r);
		

 		/**
		 * @brief       @see ReconStrategy::Name(const char*)
		 *
		 * @param name  Name
		 */
		void Name (const char* name);
		
		
 		/**
		 * @brief       @see ReconStrategy::Name()
		 *
		 * @return      Name
		 */
		const char* Name ();


	private:

		
		ReconStrategy*  m_strategy;   /**< Active strategy           */
		void*           m_dlib;       /**< Handle on startegy module */
		

		/**
		 * @brief        Default Constructor.
		 */
		ReconContext     ();
		
		
	};
	
}
#endif 
