#include "NuFFT.h"

NuFFT::NuFFT () {

	Attribute("dim", &m_dim);

	//N[0] = m_raw->dims[0];
	//N[1] = m_raw->dims[1];

	Attribute("M", &m_M);
	Attribute("b", &m_b);
	Attribute("q", &m_q);


}

RRSModule::error_code
NuFFT::Process () {

	RRSModule::error_code error = OK;

	return error;

};

