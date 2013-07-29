/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
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
#include "Algos.hpp"
#include "Creators.hpp"

using namespace RRStrategy;

#define clear(X) Free("X")

std::string sides[3] = {"Nx", "Ny", "Nz"};


NuFFT::NuFFT () : nfft(0) {

}


NuFFT::~NuFFT () {

	this->Finalise();

}


error_code 
NuFFT::Finalise () {

	return OK;

}



error_code
NuFFT::Init () {

	error_code error = OK; 
	m_initialised               = false;

	int shots, M, dim, maxit, m, N[4];
	double epsilon, alpha;
 

	for (int i = 0; i < 4; i++)
		N[i] = 1; 

	// Dimensions ---------------------------

	Attribute("dim",       &dim);

	for (int i = 0; i < dim; i++)
		Attribute (sides[i].c_str(),       &N[i]);

	Attribute("M",         &M);
	Attribute("shots",     &shots);

	// --------------------------------------

	Attribute("maxit",   &maxit);
	Attribute("epsilon", &epsilon);

	// Oversampling -------------------------

	m           = 1;
	alpha       = 1.0;

	Attribute ("m",       &m);
	Attribute ("alpha",   &alpha);

	// --------------------------------------

	Matrix<size_t> ms (dim,1);
	for (size_t i = 0; i < (size_t)dim; i++)
		ms[i] = N[i];

	nfft = new NFFT<float> (ms, M * shots, m, alpha);

	AddMatrix<cxdb> ("img");

	m_initialised = true;

	return error;

}

error_code 
NuFFT::Prepare () {

	error_code error = OK;

	nfft->KSpace  (Get<double> ("kspace"));
	nfft->Weights (Get<double> ("weights"));

	clear (kspace);
	clear (weights);

	return error;

}


error_code
NuFFT::Process () {

    SimpleTimer st ("SENSE");

    Matrix<cxdb>& out = Get<cxdb> ("img");
    Matrix<cxdb>& in = Get<cxdb> ("data");

    out = *nfft ->* in;

	clear (data);
	
	return OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new NuFFT;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

