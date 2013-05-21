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

#ifndef __TOKENIZER_HPP__
#define __TOKENIZER_HPP__

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <assert.h>

static inline std::vector<std::string>
Split     (const std::string& str, const std::string& dlm) {

	assert (dlm.size() > 0);

	std::vector<std::string> sv;
	size_t  start = 0, end = 0;
	
	while (end != std::string::npos) {
		
		end = str.find (dlm, start);
		
		// If at end, use length=maxLength.  Else use length=end-start.
		sv.push_back(str.substr(start, (end == std::string::npos) ? std::string::npos : end - start));
		
		// If at end, use start=maxSize.  Else use start=end+delimiter.
		start = ((end > (std::string::npos - dlm.size())) ? std::string::npos : end + dlm.size());

	}

	return sv;

}

#endif 
