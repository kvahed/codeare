#include "SystemCmd.hpp"

using namespace RRStrategy;

RRSModule::error_code 
SystemCmd::Process () {

	std::stringstream invoce;
	invoce << Attribute ("cmd") << " " << Attribute ("args") << " " << Attribute ("uid");

	if (std::system(NULL)) 
		puts ("Ok");
	else 
		return FATAL_SYSTEM_CALL;

	return (RRSModule::error_code ) system (invoce.str().c_str());

}
