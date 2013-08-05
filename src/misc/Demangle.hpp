#ifndef __DEMANGLE_HPP__
#define __DEMANGLE_HPP__

#include "../config.h"

#ifdef HAVE_DEMANGLE
#include <cxxabi.h>

static inline const std::string
demangle(const char* name) {
    int status = -4;
    char* res = abi::__cxa_demangle(name, NULL, NULL, &status);
    const char* const demangled_name = (status==0)?res:name;
    string ret_val(demangled_name);
    free(res);
    return ret_val;
}

#else

static inline const std::string
demangle (const char* symbol) {
    return std::string(symbol); 
}

#endif



#endif //__DEMANGLE_HPP__
