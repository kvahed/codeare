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

#ifndef __MEDIAN_FILTER_OMP_HPP__
#define __MEDIAN_FILTER_OMP_HPP__

#include "ReconStrategy.hpp"

using namespace RRServer;

namespace RRStrategy {

	/**
	 * @brief Median filter with OpenMP support
	 */
	class MedianFilter_OMP : public ReconStrategy {
		
		
	public:
		
		/**
		 * @brief Default constructor
		 */
		MedianFilter_OMP  () {};
		
		/**
		 * @brief Default destructor
		 */
		virtual 
		~MedianFilter_OMP () {};
		
		/**
		 * @brief Apply Median filter to image space
		 */
		virtual RRSModule::error_code
		Process ();
		
	};

}
#endif /* __MEDIAN_FILTER_OMP_HPP__ */
