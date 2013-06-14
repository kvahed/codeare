#ifndef __DEMANGLE_HPP__
#define __DEMANGLE_HPP__

#include "../config.h"

#ifdef HAVE_CXXABI_H
#include <cxxabi.h>

static inline std::string
demangle (const char* symbol) {

	size_t size;
	int    status;
	char   temp[128];
	char*  demangled;

	//first, try to demangle a c++ name
	if (1 == sscanf(symbol, "%*[^(]%*[^_]%127[^)+]", temp))
		if (NULL != (demangled = abi::__cxa_demangle(temp, NULL, &size, &status))) {
			std::string result(demangled);
			free(demangled);
			return result;
		}
    
	//if that didn't work, try to get a regular c symbol
	if (1 == sscanf(symbol, "%127s", temp))
		return temp;

}

#else

static inline std::string
demangle (const char* symbol) {
    return std::string(symbol); 
}

#endif



#endif //__DEMANGLE_HPP__
