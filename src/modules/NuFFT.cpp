#include "NuFFT.h"

RRSModule::error_code
NuFFT::Process () {

	m_dim  = 2;

	for (int i = 0; i < 2; i++) {
		N [i] = 0;
		Nk[i] = 0;
	}

	m_M    = 2;
	m_b    = 2.4;
	m_q    = 32;


};

