#include "Toolbox.hpp"
#include <assert.h>


Toolbox* Toolbox::m_instance = 0;


Toolbox* 
Toolbox::Instance() {

    if (m_instance == 0)
        m_instance = new Toolbox();

	return m_instance;

}


Toolbox::~Toolbox () {
     
	m_instance=0;

}


void 
Toolbox::Split (std::vector<std::string>& sv, const std::string& str, const std::string& dlm) const {
	
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

	
