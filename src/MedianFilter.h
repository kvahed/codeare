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

#ifndef __MEDIAN_FILTER_H__
#define __MEDIAN_FILTER_H__

#include "ReconStrategy.h"

/**
 * @brief Median filter
 */
class MedianFilter : public ReconStrategy {


public:

	/**
	 * @brief Default constructor
	 */
	MedianFilter  () {};
	
	/**
	 * @brief Default destructor
	 */
	virtual 
	~MedianFilter () {};

	/**
	 * @brief Apply Median filter to image space
	 */
	virtual RRSModule::error_code
	ProcessData ();


private:	

	int m_ww; /**> @brief filter window width*/
	int m_wh; /**> @brief filter window height*/

};

#endif /* __MEDIAN_FILTER_H__ */
