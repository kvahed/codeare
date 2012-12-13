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
#include "IO.hpp"

using namespace RRStrategy;

std::string sides[3] = {"Nx", "Ny", "Nz"};


NuFFT::NuFFT () {

}


NuFFT::~NuFFT () {

	this->Finalise();

}


RRSModule::error_code 
NuFFT::Finalise () {

	//if (m_initialised)
	//	delete nfft;

	return RRSModule::OK;

}



RRSModule::error_code
NuFFT::Init () {

	RRSModule::error_code error = OK; 
	m_initialised               = false;

	int shots, M, dim, maxit, m, N[4], n[4];
	double epsilon, alpha;
 

	for (int i = 0; i < 4; i++) {
		N[i] = 1; 
		n[i] = 1;
	}

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

	for (int i = 0; i < dim; i++)
		n[i] = ceil (N[i]*alpha);

	// --------------------------------------

	Matrix<size_t> ms (dim,1);
	for (size_t i = 0; i < dim; i++)
		ms[i] = N[i];

	nfft = new NFFT<float> (ms, M * shots, m, alpha);

	Matrix<cxdb>& img = AddMatrix 
		("img", (Ptr<Matrix<cxdb> >) NEW (Matrix<cxdb> (N[0],N[1],N[2])));

	m_initialised = true;

	return error;

}

RRSModule::error_code 
NuFFT::Prepare () {

	RRSModule::error_code error = OK;

	nfft->KSpace  (Get<double> ("kspace"));
	nfft->Weights (Get<double> ("weights"));

	Free<double> ("kspace");
	Free<double> ("weights");

	return error;

}


RRSModule::error_code
NuFFT::Process () {

	printf ("Processing NuFFT ...\n");
	ticks start = getticks();
	
	Get<cxdb> ("img") = nfft->Adjoint(Get<cxdb> ("data"));
	
	printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), start)/Toolbox::Instance()->ClockRate());

	return OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new NuFFT;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

