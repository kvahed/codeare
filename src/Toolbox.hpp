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

#ifndef __TOOLBOX_HPP__
#define __TOOLBOX_HPP__

#include <boost/any.hpp>

#include <vector>
#include <string>

using namespace std;

/**
 * @brief  A toolbox for some static stuff
 */
class Toolbox {

public:

	/**
	 * @brief       Singleton destructor
	 */
	~Toolbox();


    /**
     * @brief       Singleton instance
     */
    static Toolbox*  
    Instance        ();


	/**
	 * @brief       CPU clock rate.
	 *              Humble abuse of FFTW cycle for timing information.
	 *
	 * @return      clock rate
	 */
	double 
	ClockRate       () const ;	

	boost::any   void_any;
	std::string  void_str;
	
		
private:


	/**
	 * @brief      Hide constructor for singletons
	 */
	Toolbox   () {};


	static Toolbox*    m_instance;           /**< @brief Single instance */

};

#endif
