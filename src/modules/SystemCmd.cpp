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

#include "SystemCmd.hpp"

using namespace RRStrategy;

codeare::error_code 
SystemCmd::Process () {

	std::stringstream invoce;
	invoce << Attribute ("cmd") << " " << Attribute ("args") << " " << Attribute ("uid");

	if (std::system(NULL)) 
		puts ("Ok");
	else 
	  return codeare::FATAL_SYSTEM_CALL;

	return (codeare::error_code ) system (invoce.str().c_str());

}
