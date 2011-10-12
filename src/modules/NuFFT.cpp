/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301  USA
 */

#include "NuFFT.hpp"

using namespace RRStrategy;

std::string sides[3] = {"Nx", "Ny", "Nz"};


NuFFT::NuFFT () {

}


NuFFT::~NuFFT () {

	this->Finalise();

}


RRSModule::error_code 
NuFFT::Finalise () {

	if (m_initialised)
	  nfft::finalize (&m_fplan, &m_iplan);

	delete [] m_N;
	delete [] m_n;

	return RRSModule::OK;

}



RRSModule::error_code
NuFFT::Init () {

	RRSModule::error_code error = OK; 
	m_initialised               = false;

	m_N   = new int[3];
	m_n   = new int[3];

	for (int i = 0; i < 3; i++) {
		m_N[i] = 0; 
		m_n[i] = 0;
	}

	// Dimensions ---------------------------

	Attribute("dim",       &m_dim);

	for (int i = 0; i < m_dim; i++)
		Attribute (sides[i].c_str(),       &m_N[i]);

	Attribute("M",         &m_M);
	Attribute("shots",     &m_shots);

	// --------------------------------------

	Attribute("maxit",   &m_maxit);
	Attribute("epsilon", &m_epsilon);

	// Oversampling -------------------------

	int      m           = 1;
	double   alpha       = 1.0;

	Attribute ("m",       &m);
	Attribute ("alpha",   &alpha);

	for (int i = 0; i < m_dim; i++)
		m_n[i] = ceil (m_N[i]*alpha);

	// --------------------------------------

	printf ("  intialising nfft::init (%i, {%i, %i, %i}, %i, {%i, %i, %i}, %i, *, *, %.9f)\n", 
			m_dim, 
			m_N[0], m_N[1], m_N[2],
			m_M * m_shots,
			m_n[0], m_n[1], m_n[2],
			m,
			m_epsilon);

	nfft::init (m_dim, m_N, m_M*m_shots, m_n, m, &m_fplan, &m_iplan, m_epsilon);

	m_initialised = true;

	return error;

}

RRSModule::error_code
NuFFT::Process () {

	// Some variables
	RRSModule::error_code error = OK;

	Matrix<cplx>&   data    = GetCplx("data");
	Matrix<double>& kspace  = GetReal("kspace");
	Matrix<double>& weights = GetReal("weights");

	printf ("Processing NuFFT ...\n");
	ticks start = getticks();

	// Copy data from incoming matrix to the nufft input array
	for (int i = 0; i < data.Size(); i++) {
		m_iplan.y[i][0] = data[i].real();
		m_iplan.y[i][1] = data[i].imag();
	}

	// Copy k-space and weights to allocated memory
	memcpy (&(m_fplan.x[0]),  &kspace[0],  kspace.Size()*sizeof(double));
	memcpy (&(m_iplan.w[0]), &weights[0], weights.Size()*sizeof(double));

	// Assign weights and precompute PSI
	nfft::weights (&m_fplan, &m_iplan);
	nfft::psi (&m_fplan);

	// Inverse FT
	nfft::ift (&m_fplan, &m_iplan, m_maxit, m_epsilon);

	// Resize data for output
	for (int i = 0; i < INVALID_DIM; i++)
		data.Dim(i) = 1;
	for (int i = 0; i < m_dim; i++)
		data.Dim(i) = m_N[i];
	data.Reset();

	// Copy back reconstructed image to outgoing matrix
	for (int i = 0; i < data.Size(); i++)
		data[i] = cplx (m_iplan.f_hat_iter[i][0], m_iplan.f_hat_iter[i][1]);
	
	printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), start) / Toolbox::Instance()->ClockRate());
	
	return error;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new NuFFT;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

