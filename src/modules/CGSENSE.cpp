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

#include "CGSENSE.hpp"
#include "Noise.hpp"
#include "Toolbox.hpp"
#include "Lapack.hpp"
#include "Creators.hpp"

#include <math.h>


#include <vector>

std::string sides[3] = {"nx", "ny", "nz"};

using namespace RRStrategy;


CGSENSE::~CGSENSE () {

	this->Finalise();

}


error_code 
CGSENSE::Finalise () {

	if (m_initialised)
		delete m_ncs;

	return OK;

}


error_code 
CGSENSE::Init() {

	printf ("Intialising CG-SENSE ...\n");

	error_code error = OK; 

	m_initialised = false;

	// Some defaults ------------------------
	m_testcase = 0;
	m_verbose  = 0;
	m_noise    = 0;
	m_lambda   = 5.0e-2;

	// --------------------------------------

	// Tikhonov ----------------------------

	Attribute ("lambda",  &m_lambda);
	// --------------------------------------

	// Verbosity ----------------------------

	Attribute ("testcase",  &m_testcase);
	printf ("  test case: %i \n", m_testcase);
	// --------------------------------------

	// Verbosity ----------------------------

	Attribute ("verbose",   &m_verbose);
	printf ("  verbose feedback: %i \n", m_verbose);
	// --------------------------------------

	// Noise --------------------------------

	Attribute ("noise",   &m_noise);
	printf ("  gaussian white noise (normalised): %.9f \n", m_noise);
	// --------------------------------------

	// CG convergence and break criteria ----

	Attribute ("cgeps",   &m_cgeps);
	Attribute ("cgmaxit", &m_cgmaxit);
	printf ("  maximum #iterations: %i \n", m_cgmaxit);
	printf ("  convergence criterium: %.9f \n", m_cgeps);
	// --------------------------------------

	// iNFFT convergence and break criteria -
	
	Attribute ("ftmaxit", &m_ftmaxit);
	Attribute ("fteps",   &m_fteps);
	// --------------------------------------

	printf ("... done.\n\n");

	return error;

}


error_code
CGSENSE::Prepare () {

	error_code error = OK;

	Matrix<float>& weights = Get<float>("weights");
	Matrix<cxfl>& sens = Get<cxfl>("sens");
	Matrix<float>& kspace = Get<float>("kspace");

	size_t nk = numel(weights);

	Params cgp;
	cgp.Set("sens_maps", std::string("sens"));
	cgp.Set("sens_maps", std::string("weights"));
	cgp.Set("sens_maps", std::string("kspace"));


	m_ncs = new NCSENSE<float>
		(sens, nk, m_cgeps, m_cgmaxit, m_lambda, m_fteps, m_ftmaxit);

	size_t dim = ndims(sens) - 1;

	// Outgoing images
	Matrix<cxfl>& image = AddMatrix 
		("image", (Ptr<Matrix<cxfl> >) NEW (Matrix<cxfl>(size(sens,0), size(sens,1), (dim == 3) ? size(sens,2) : 1)));

	m_ncs->KSpace (kspace);
	m_ncs->Weights (weights);
	
	Free ("weights");
	Free ("kspace");
	
	m_initialised = true;

	return error;

}


error_code
CGSENSE::Process () {

	error_code error = OK;

	ticks cgstart = getticks();
	
	printf ("Processing CGSENSE ...\n");

	Get<cxfl>("image") = m_ncs->Adjoint (Get<cxfl>("data"), Get<cxfl>("sens"));

	Free ("data");

	printf ("... done. WTime: %.4f seconds.\n\n", elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());

	return error;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}





