#include "LocalConnector.hpp"

namespace RRClient {

LocalConnector::~LocalConnector ()               {
	
	FunctorContainer::Finalise ();
	
}


error_code
LocalConnector::CleanUp () {
	
	FunctorContainer::Finalise ();
	
}

error_code
LocalConnector::Init (const char* name) {

	std::stringstream  temp;
	temp << GetConfig();
	FunctorContainer::config  (temp.str().c_str());
		
	return FunctorContainer::Init (name);

}


error_code
LocalConnector::Finalise (const char* name) {
	
	return FunctorContainer::Finalise (name);
	
}


error_code
LocalConnector::Process  (const char* name)       {
	
	return FunctorContainer::Process (name);
	
}


error_code
LocalConnector::Prepare  (const char* name)       {
	
	return FunctorContainer::Prepare (name);
	
}


}
