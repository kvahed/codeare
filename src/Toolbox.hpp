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

#ifndef __TOOLBOX_HPP__
#define __TOOLBOX_HPP__

#include <vector>
#include <string>

using namespace std;

static std::string bars   = "***************************************************";
static std::string blancs = "                                                   ";
static std::string spear  = ">";

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
	 * @brief       String splitter
	 *
	 * @param  sv   Splitted parts
	 * @param  str  String for splitting
	 * @param dlm   Delimiter
	 */
	void 
	Split           (std::vector<std::string>& sv, const std::string& str, const std::string& dlm) const;


	/**
	 * @brief       CPU clock rate.
	 *              Humble abuse of FFTW cycle for timing information.
	 *
	 * @return      clock rate
	 */
	double 
	ClockRate       () const ;	


	
	/**
	 * @brief       Command line progress bar to cout
	 *
	 * @param  pre  Text before bar
	 * @param  post Test after bar
	 * @param  p    Percent
	 */ 
	void 
	ProgressBar    (const std::string& pre, const std::string& post, const short& p) const;

		
private:


	/**
	 * @brief      Hide constructor for singletons
	 */
	Toolbox   () {};

	static Toolbox*    m_instance;           /**< @brief Single instance */

};

#endif
