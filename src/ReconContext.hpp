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

#include <vector>

#include "ReconStrategy.hpp"

namespace RRServer {

	/**
	 * @brief Context of a reconstruction method
	 */
	class ReconContext {
		
		
		
	public:
		
		
		/**
		 * @brief Default Constructor
		 */
		ReconContext () {}
		
		
		/**
		 * @brief Invoce destruction on my startegy and exit
		 */ 
		~ReconContext ();
		
		
		/**
		 * @brief Construct with a strategy
		 */
		ReconContext (const char* name);
		
		
		/**
		 * @brief get active startegy
		 */
		inline ReconStrategy*
			Strategy     () {
			return m_strategy;
		}
		
		/**
		 * @brief Process data with given strategy
		 */
		RRSModule::error_code
			Process () {
			return m_strategy->Process();
		}
		
		
		/**
		 * @brief Process data with given strategy
		 */
		RRSModule::error_code
			Init () {
			return m_strategy->Init();
		}
		
		
		/**
		 * @brief Process data with given strategy
		 */
		RRSModule::error_code
			Finalise () {
			return m_strategy->Finalise();
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		SetConfig (const char* cstr) {
			m_strategy->SetConfig(cstr);
		}
		
		
		/**
		 * @brief Process data with given strategy
		 */
		void
		ReadConfig (const char* fname) {
			m_strategy->ReadConfig(fname);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		GetConfig (char* cstr) {
			m_strategy->GetConfig(cstr);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		SetRaw (const raw_data* r) {
			m_strategy->SetRaw(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		SetRaw (const Matrix<raw>* r) {
			m_strategy->SetRaw(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		GetRaw (Matrix<raw>* r) {
			m_strategy->GetRaw(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		GetRaw (raw_data* r) {
			m_strategy->GetRaw(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		SetRHelper (const raw_data* r) {
			m_strategy->SetRHelper(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		SetRHelper (const Matrix<raw>* r) {
			m_strategy->SetRHelper(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		GetRHelper (Matrix<raw>* r) {
			m_strategy->GetRHelper(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		GetRHelper (raw_data* r) {
			m_strategy->GetRHelper(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		SetHelper (const helper_data* r) {
			m_strategy->SetHelper(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		SetHelper (const Matrix<double>* r) {
			m_strategy->SetHelper(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		GetHelper (Matrix<double>* r) {
			m_strategy->GetHelper(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		GetHelper (helper_data* r) {
			m_strategy->GetHelper(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		SetKSpace (const helper_data* r) {
			m_strategy->SetKSpace(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		SetKSpace (const Matrix<double>* r) {
			m_strategy->SetKSpace(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		GetKSpace (Matrix<double>* r) {
			m_strategy->GetKSpace(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		GetKSpace (helper_data* r) {
			m_strategy->GetKSpace(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		SetPixel (const pixel_data* r) {
			m_strategy->SetPixel(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		SetPixel (const Matrix<short>* r) {
			m_strategy->SetPixel(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		GetPixel (Matrix<short>* r) {
			m_strategy->GetPixel(r);
		}
		

		/**
		 * @brief Process data with given strategy
		 */
		void
		GetPixel (pixel_data* r) {
			m_strategy->GetPixel(r);
		}
		

		void Name (const char* name) { m_strategy->Name(name);}
		
		const char* Name () {return m_strategy->Name();}
		
	private:
		
		ReconStrategy*            m_strategy;   /**< Active strategy           */
		void*                     m_dlib;       /**< Handle on startegy module */
		
	};
	
}
#endif 
