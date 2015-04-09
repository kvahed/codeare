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
#include "Access.hpp"

using namespace RRStrategy;

#define clear(X) Free("X")

std::string sides[3] = {"Nx", "Ny", "Nz"};


NuFFT::NuFFT () : m_dft_3rd_dim(false), m_test_case(false) {}


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
	float epsilon, alpha;
 

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
	Attribute("dft_3rd_dim", &m_dft_3rd_dim);
	Attribute("test_case", &m_test_case);

	// Oversampling -------------------------

	m           = 1;
	alpha       = 1.0;

	Attribute ("m",       &m);
	Attribute ("alpha",   &alpha);

	// --------------------------------------

	Vector<size_t> ms (dim,1);
	for (size_t i = 0; i < (size_t)dim; i++)
		ms[i] = N[i];


	Matrix<float> b0 = (Exists<float>("b0")==codeare::OK) ?
			Get<float>("b0") : Matrix<float>(1);
	Matrix<float> timing = (Exists<float>("timing")==codeare::OK) ?
			Get<float>("timing") : Matrix<float>(1);

	Params p;
	p["alpha"]  = alpha;
	p["m"]      = m;
	p["b0"]     = b0;
	p["timing"] = timing;
	p["nk"]     = M*shots;
	p["3rd_dim_cart"] = m_dft_3rd_dim;
	p["imsz"]   = ms;

	ft = NFFT<cxfl> (p);

	m_initialised = true;

	return error;

}

codeare::error_code 
NuFFT::Prepare () {

	codeare::error_code error = codeare::OK;

	ft.KSpace  (Get<float> ("kspace"));
	ft.Weights (Get<float> ("weights"));

	std::cout << ft << std::endl;

	clear (kspace);
	clear (weights);

	return error;

}


codeare::error_code
NuFFT::Process () {

	Matrix<cxfl> img;
	const Matrix<cxfl>& data = (m_test_case) ?
			ft * phantom<cxfl>(ft.ImageSize()) : Get<cxfl> ("data");
	img = ft ->* data;

    Add ("img", img);
	return codeare::OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new NuFFT;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

