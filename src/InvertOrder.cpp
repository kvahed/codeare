#include "InvertOrder.h"

RRSModule::error_code
InvertOrder::ProcessData () {

	/*	RRSModule::floats tmpabs (m_raw.dabs);
	RRSModule::floats tmparg (m_raw.darg);
	
	int     len = tmpabs.length();
	
	for (int i = 0; i < len; i++) {
		tmpabs[i] = m_raw.dabs [len-1-i];
		tmparg[i] = m_raw.darg [len-1-i];
	}
	
	m_raw.dabs = tmpabs;
	m_raw.darg = tmparg;*/
	
	return RRSModule::OK;

}
