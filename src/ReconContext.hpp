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

#ifndef __RECONCONTEXT_HPP__
#define __RECONCONTEXT_HPP__

#include "ReconStrategy.hpp"
#include <vector>

namespace RRServer {


	/**
	 * @brief Context of a reconstruction method. Abstraction layer to algorithm backends. 
	 */
	class ReconContext {
		
		
		
	public:
		
		
		/**
		 * @brief        Default Constructor.
		 */
		ReconContext     ();
		
		
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
		 * @brief        @see ReconStrategy::SetRaw(const raw_data*)
		 *
		 * @param  r     Complex data sequence
		 */
		void
		SetRaw          (const raw_data* r);
		

		/**
		 * @brief       @see ReconStrategy::SetRaw(const Matrix<raw>*)
		 *
		 * @param  r    Complex data matrix
		 */
		void
		SetRaw          (const Matrix<raw>* r);
		

		/**
		 * @brief       @see ReconStrategy::GetRaw(raw_data*)
		 *
		 * @param  r    Complex data sequence
		 */
		void
		GetRaw          (raw_data* r);
		

		/**
		 * @brief       @see ReconStrategy::GetRaw(Matrix<raw>*)
		 *
		 * @param  r    Complex data matrix
		 */
		void
		GetRaw          (Matrix<raw>* r);
		

		/**
		 * @brief       @see ReconStrategy::SetRHelper(const raw_data*)
		 *
		 * @param  r    Complex data sequence
		 */
		void
		SetRHelper (const raw_data* r);
		

		/**
		 * @brief       @see ReconStrategy::SetRHelper(const Matrix<raw>*)
		 *
		 * @param  r    Complex data matrix
		 */
		void
		SetRHelper (const Matrix<raw>* r);
		

		/**
		 * @brief       @see ReconStrategy::GetRHelper(raw_data*)
		 *
		 * @param  r    Complex data sequence
		 */
		void
		GetRHelper (raw_data* r);
		

		/**
		 * @brief       @see ReconStrategy::GetRHelper(Matrix<raw>*)
		 *
		 * @param  r    Complex data matrix
		 */
		void
		GetRHelper (Matrix<raw>* r);
		

		/**
		 * @brief       @see ReconStrategy::SetHelper(const helper_data*)
		 *
		 * @param  r    Real data sequence
		 */
		void
		SetHelper (const helper_data* r);
		

		/**
		 * @brief       @see ReconStrategy::SetHelper(const Matrix<double>*)
		 *
		 * @param  r    Real data matrix
		 */
		void
		SetHelper (const Matrix<double>* r);
		

		/**
		 * @brief       @see ReconStrategy::SetHelper(const Matrix<double>*)
		 *
		 * @param  r    Real data matrix
		 */
		void
		GetHelper (Matrix<double>* r);
		

		/**
		 * @brief       @see ReconStrategy::GetHelper(helper_data*)
		 *
		 * @param  r    Real data sequence
		 */
		void
		GetHelper (helper_data* r);
		

		/**
		 * @brief       @see ReconStrategy::SetKSpace(helper_data*)
		 *
		 * @param  r    Real data sequence
		 */
		void
		SetKSpace (const helper_data* r);
		

		/**
		 * @brief       @see ReconStrategy::SetKSpace(const Matrix<double>*)
		 *
		 * @param  r    Real data matrix
		 */
		void
		SetKSpace (const Matrix<double>* r);
		

		/**
		 * @brief       @see ReconStrategy::GetKSPace(Matrix<double>*)
		 *
		 * @param  r    Real data sequence
		 */
		void
		GetKSpace (Matrix<double>* r);
		

		/**
		 * @brief       @see ReconStrategy::GetKSpace(helper_data*)
		 *
		 * @param  r    Real data sequence
		 */
		void
		GetKSpace (helper_data* r);
		

		/**
		 * @brief       @see ReconStrategy::SetPixel(const pixel_data*)
		 *
		 * @param  r    Pixel data sequence
		 */
		void
		SetPixel (const pixel_data* r);
		

		/**
		 * @brief       @see ReconStrategy::SetPixel(const Matrix<short>*)
		 *
		 * @param  r    Pixel data matrix
		 */
		void
		SetPixel (const Matrix<short>* r);
		

		/**
		 * @brief       @see ReconStrategy::GetPixel(Matrix<short>*)
		 *
		 * @param  r    Pixel data matrix
		 */
		void
		GetPixel (Matrix<short>* r);
		

		/**
		 * @brief       @see ReconStrategy::GetPixel(pixel_data*)
		 *
		 * @param  r    Pixel data sequence
		 */
		void
		GetPixel (pixel_data* r);
		

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
		
	};
	
}
#endif 
