#include "NuFFT.hpp"

using namespace RRStrategy;

NuFFT::NuFFT () {

	Attribute("dim",     &m_dim);

	Attribute("N",       &m_N);
	Attribute("M",       &m_M);
	Attribute("iter",    &m_iter;)
	Attribute("epsilon", &m_epsilon;)

}

RRSModule::error_code
NuFFT::Process () {

	RRSModule::error_code error = OK;

	Matrix<raw> in = m_raw;

	for (int i = 0; i < m_dim; i++)
		m_raw.Dim(i) = m_N;

	for (int i = m_dim; i < INVALID_DIM; i++)
		m_raw.Dim(i) =1 ;
	
	m_raw.Reset();

	double*     ftin     = new double[2 *  in->Size()]; 

	for (int i = 0; i < in->Size(); i++) {
		ftin[2*i  ] = (in->at(i)).real();
		ftin[2*i+1] = (in->at(i)).imag();
	}
		
	nfft::ift (m_fplan, m_iplan, ftin, ftout, m_iter, m_epsilon);

	m_raw->at(i) = raw(ftout[2*i] * sens.real(), ftout[2*i+1] * -sens.imag()); 

	delete [] ftout;
	delete [] ftin;

	return error;

}

