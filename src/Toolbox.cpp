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

#include "Toolbox.hpp"
#include "config.h"

#include <assert.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <string>

std::string exec(char* cmd) {

  FILE* pipe = popen(cmd, "r");
  if (!pipe) return "ERROR";
  char buffer[128];
  std::string result = "";

  while(!feof(pipe))
    if(fgets(buffer, 128, pipe) != NULL)
      result += buffer;

  pclose(pipe);
  return result;

}

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
	std::string mhzstr = exec("lscpu | grep \"CPU MHz\"|awk '{print $3}'");
	float mhz = atof(mhzstr.c_str());
	return 1000000.0 * mhz;

#endif
}

	
