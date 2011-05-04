#include "NuFFT.hpp"

using namespace RRStrategy;

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

	Matrix<raw> tmp = m_raw;

	/*	// Clear outgoing container
	out->Reset();

	// Some dimensions
	int         ncoils   = sm->Dim(CHA);
	int         nsamples = in->Size(); 
	int         ndim     = out->Dim(COL);

	// Container for FT input
	double*     ftin     = new double[2 *  in->Size()/ncoils]; 

	// Loop over coils, Inverse FT every signal in *in
	// and Sum of squares of images
	for (int j = 0; j < ncoils; j++) {
		
		double* ftout = new double[2 * out->Size()];
		int     pos   = j*in->Dim(COL)*in->Dim(LIN);

		for (int i = 0; i < in->Size()/ncoils; i++) {
			ftin[2*i  ] = (in->at(pos + i)).real();
			ftin[2*i+1] = (in->at(pos + i)).imag();
		}
		
		nfft::ift (np, spc, ftin, ftout, maxit, epsilon);

		for (int i = 0; i < out->Size(); i++) {
			raw sens = sm->at(pos + i);
			out->at(i) += raw(ftout[2*i] * sens.real(), ftout[2*i+1] * -sens.imag());
		}

		delete [] ftout;

	}
	
	delete [] ftin;
	*/
	return error;

}

