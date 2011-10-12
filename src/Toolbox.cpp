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

#include "Toolbox.hpp"
#include "config.h"

#include <assert.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>

Toolbox* Toolbox::m_instance = 0;


Toolbox::~Toolbox  () {
     
	m_instance=0;

}


Toolbox* 
Toolbox::Instance  () {

    if (m_instance == 0)
        m_instance = new Toolbox();

	return m_instance;

}


void 
Toolbox::Split     (std::vector<std::string>& sv, const std::string& str, const std::string& dlm) const {
	
	assert (dlm.size() > 0);
	
	size_t  start = 0, end = 0;
	
	while (end != std::string::npos) {
		
		end = str.find (dlm, start);
		
		// If at end, use length=maxLength.  Else use length=end-start.
		sv.push_back(str.substr(start, (end == std::string::npos) ? std::string::npos : end - start));
		
		// If at end, use start=maxSize.  Else use start=end+delimiter.
		start = ((end > (std::string::npos - dlm.size())) ? std::string::npos : end + dlm.size());

	}

}



double 
Toolbox::ClockRate () const {
	
#if defined(HAVE_MACH_ABSOLUTE_TIME)
	
	// OSX
	
	uint64_t freq = 0;
	size_t   size = sizeof(freq);
	
	if (sysctlbyname("hw.tbfrequency", &freq, &size, NULL, 0) < 0)
		perror("sysctl");
	
	return freq;
	
#else
	
	// LINUX
	
	FILE* scf;
	std::string fname = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq";
	int   freq = 1;
	
	scf = fopen(fname.c_str(), "rb");
	
	if (scf != NULL) {
		int read = fscanf(scf,"%i",&freq);
#ifdef VERBOSE
		printf ("Read %i from %s\n", read, fname.c_str());
#endif
		fclose(scf);
	}
	
	return 1000.0 * freq;
	
#endif
}

	
inline void 
Toolbox::ProgressBar (const std::string& pre, const std::string& post, const short& p) const {
	
	assert (p >=  0);
	assert (p <=100);
	
	std::cout << "\r";
	std::cout << pre.c_str();
	std::cout << " | "; 
	std::cout << bars.substr(0, p/2) << "> " <<  blancs.substr(0, 50-p/2) << "| " << std::setw(3) << std::setfill(' ') << p << "% done";
	
}
