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


NuFFT::NuFFT () {

}


NuFFT::~NuFFT () {

	this->Finalise();

}


codeare::error_code 
NuFFT::Finalise () {

	return codeare::OK;

}



codeare::error_code
NuFFT::Init () {

	codeare::error_code error = codeare::OK; 
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

	ft = NFFT<double> (ms, M * shots, m, alpha);

	m_initialised = true;

	return error;

}

codeare::error_code 
NuFFT::Prepare () {

	codeare::error_code error = codeare::OK;

	ft.KSpace  (Get<double> ("kspace"));
	ft.Weights (Get<double> ("weights"));

	clear (kspace);
	clear (weights);

	return error;

}


codeare::error_code
NuFFT::Process () {

    SimpleTimer st ("NuFFT");

    Matrix<cxdb> img =  ft ->* Get<cxdb> ("data");
    wspace.Add ("img", img);

	clear (data);
	
	return codeare::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new NuFFT;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

