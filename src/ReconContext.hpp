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
		 * @brief        @see ReconStrategy::SetCplx(const std::string name, const cplx_data*)
		 *
		 * @param  name  Name
		 * @param  r     Complex data sequence
		 */
		void
		SetCplx          (const std::string name, const cplx_data& r);
		

		/**
		 * @brief       @see ReconStrategy::SetCplx(const std::string name, const Matrix<cplx>*)
		 *
		 * @param  name  Name
		 * @param  r    Complex data matrix
		 */
		void
		SetCplx          (const std::string name, Matrix<cplx>& r);
		

		/**
		 * @brief       @see ReconStrategy::GetCplx(const std::string name, cplx_data&)
		 *
		 * @param  name  Name
		 * @param  r    Complex data sequence
		 */
		void
		GetCplx          (const std::string name, cplx_data& r);
		

		/**
		 * @brief       @see ReconStrategy::GetCplx(const std::string name, Matrix<cplx>*)
		 *
		 * @param  name  Name
		 * @param  r    Complex data matrix
		 */
		void
		GetCplx          (const std::string name, Matrix<cplx>& r);
		

		/**
		 * @brief       @see ReconStrategy::SetReal(const std::string name, const real_data*)
		 *
		 * @param  name  Name
		 * @param  r    Real data sequence
		 */
		void
		SetReal (const std::string name, const real_data& r);
		

		/**
		 * @brief       @see ReconStrategy::SetReal(const std::string name, const Matrix<double>*)
		 *
		 * @param  name  Name
		 * @param  r    Real data matrix
		 */
		void
		SetReal (const std::string name, Matrix<double>& r);
		

		/**
		 * @brief       @see ReconStrategy::SetReal(const std::string name, const Matrix<double>*)
		 *
		 * @param  name  Name
		 * @param  r    Real data matrix
		 */
		void
		GetReal (const std::string name, Matrix<double>& r);
		

		/**
		 * @brief       @see ReconStrategy::GetReal(const std::string name, real_data*)
		 *
		 * @param  name  Name
		 * @param  r    Real data sequence
		 */
		void
		GetReal (const std::string name, real_data& r);
		

		/**
		 * @brief       @see ReconStrategy::SetPixel(const std::string name, const pixel_data*)
		 *
		 * @param  name  Name
		 * @param  r    Pixel data sequence
		 */
		void
		SetPixel (const std::string name, const pixel_data& r);
		

		/**
		 * @brief       @see ReconStrategy::SetPixel(const std::string name, const Matrix<short>*)
		 *
		 * @param  name  Name
		 * @param  r    Pixel data matrix
		 */
		void
		SetPixel (const std::string name, Matrix<short>& r);
		

		/**
		 * @brief       @see ReconStrategy::GetPixel(const std::string name, Matrix<short>*)
		 *
		 * @param  name  Name
		 * @param  r    Pixel data matrix
		 */
		void
		GetPixel (const std::string name, Matrix<short>& r);
		

		/**
		 * @brief       @see ReconStrategy::GetPixel(const std::string name, pixel_data*)
		 *
		 * @param  name  Name
		 * @param  r    Pixel data sequence
		 */
		void
		GetPixel (const std::string name, pixel_data& r);
		

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
