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
#include "SimpleTimer.hpp"
#include "Algos.hpp"
#include "Creators.hpp"

#include "Print.hpp"

#include <math.h>


#include <vector>

std::string sides[3] = {"nx", "ny", "nz"};

using namespace RRStrategy;


CGSENSE::~CGSENSE () {
	this->Finalise();
}


codeare::error_code 
CGSENSE::Finalise () {
	return codeare::OK;
}


codeare::error_code 
CGSENSE::Init() {

	printf ("Intialising CG-SENSE ...\n");

	codeare::error_code error = codeare::OK; 

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

	// # threads ----------------------------

	Attribute ("threads",  &m_nthreads);
	printf ("  # of threads: %i \n", m_nthreads);
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


codeare::error_code
CGSENSE::Prepare () {

	codeare::error_code error = codeare::OK;

	Params cgp;
	cgp["sens_maps"]    = std::string("sensitivities");
    cgp["weights_name"] = std::string("weights");
    cgp["verbose"]      = m_verbose;
    cgp["ftiter"]       = (size_t) m_ftmaxit;
    cgp["cgiter"]       = (size_t) m_cgmaxit;
    cgp["cgeps"]        = m_cgeps;
    cgp["lambda"]       = m_lambda;
    cgp["np"]           = m_nthreads;

	m_ncs = NCSENSE<float>(cgp);

	m_ncs.KSpace (Get<float>("kspace"));
	m_ncs.Weights (Get<float>("weights"));
	
	Free ("weights");
	Free ("kspace");
	
	m_initialised = true;

	return error;

}


codeare::error_code
CGSENSE::Process () {

	codeare::error_code error = codeare::OK;
    const Matrix<cxfl>& sens = Get<cxfl>("sensitivities");
    Matrix<cxfl> data;

    data = (!m_testcase) ? Get<cxfl>("signals") :
        m_ncs.Trafo (phantom<cxfl>(size(sens,0)), sens);
    if (m_noise)
        data += m_noise * randn<cxfl>(size(data));

    Matrix<cxfl> img = m_ncs ->* data;
    
    wspace.Add("image", img);
	return error;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}





