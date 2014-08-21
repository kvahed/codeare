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

#ifndef __GRAPPA_HPP__
#define __GRAPPA_HPP__

#include "ReconStrategy.hpp"
#include "Algos.hpp"
#include "Toolbox.hpp"
#include "Creators.hpp"
#include "CGRAPPA.hpp"

/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

	/**
	 * @brief GRAPPA PPI reconstruction
	 */
	class GRAPPA : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		GRAPPA () : m_ft(0) {};


		/**
		 * @brief Default destructor
		 */
		virtual 
		~GRAPPA () {};
		
		
		/**
		 * @brief Do nothing 
		 */
		virtual codeare::error_code
		Process ();
		
		/**
		 * @brief Do nothing 
		 */
		virtual codeare::error_code
		Init ();
		
		/**
		 * @brief Do nothing 
		 */
		virtual codeare::error_code
		Prepare ();
		
		/**
		 * @brief Do nothing 
		 */
		virtual codeare::error_code
		Finalise ();


	private: 
		
        CGRAPPA<double>* m_ft;
        Vector<size_t> m_kernel_size, m_acceleration_factors;
        size_t m_nthreads;
        float m_lambda;
	  
	};

}

#endif /* __GRAPPA_H__ */

