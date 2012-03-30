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

#include "CGSENSE.hpp"
#include "nfftstub.h"
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


RRSModule::error_code 
CGSENSE::Finalise () {

	if (m_initialised)
		for (int i = 0; i < NTHREADS || i < m_Nc; i++)
			nfft::finalize (&m_fplan[i], &m_iplan[i]);

	return OK;

}


RRSModule::error_code 
CGSENSE::Init() {

	printf ("Intialising CG-SENSE ...\n");

	RRSModule::error_code error = OK; 

	m_initialised = false;

	// Some defaults ------------------------
	for (int i = 0; i < 3; i++) {
		m_N[i] = 1; 
		m_n[i] = 1;
	}

	m_testcase = 0;
	m_verbose  = 0;
	m_noise    = 0;
	m_dim      = 1;
	m_M        = 0;
	m_lambda   = 5.0e-2;

	// Dimensions ---------------------------

	Attribute("dim",       &m_dim);
	printf ("  dimensions: %iD \n", m_dim);


	if (m_dim < 2 || m_dim > 3) {
		printf ("%s only supports 2-3 dimensions. %i was specified.\n", Name(), m_dim);
		return UNSUPPORTED_DIMENSION;
	}

	for (int i = 0; i < m_dim; i++)
		Attribute (sides[i].c_str(),       &m_N[i]);

	for (int i = 0; i < m_dim; i++)
		if (m_N[i] < 16) {
			printf ("%s only supports image matrix sides >= 16. (%ix%ix%i) was specified.\n", Name(), m_N[0], m_N[1], m_N[2]);
			return UNSUPPORTED_IMAGE_MATRIX;
		}

	Attribute ("M",         &m_M);

	if (m_M == 0) {
		printf ("Initialising %s with %i mesurement nodes? Check configuration! FAILED: Bailing out!\n", Name(), m_M);
		return ZERO_NODES;
	}

	Attribute ("Nc",        &m_Nc);

	if (m_Nc == 0) {
		printf ("Initialising %s with %i channels? Check configuration! FAILED: Bailing out!\n", Name(), m_Nc);
		return CGSENSE_ZERO_CHANNELS;
	}

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
	
	Attribute ("maxit",   &m_maxit);
	Attribute ("epsilon", &m_epsilon);
	// --------------------------------------

	// Oversampling -------------------------

	int      m           = 1;
	double   alpha       = 1.0;

	Attribute ("m",       &m);
	Attribute ("alpha",   &alpha);

	for (int i = 0; i < m_dim; i++)
		m_n[i] = ceil (m_N[i]*alpha);
	// --------------------------------------

	// Initialise FT plans ------------------
	
	printf ("  intialising nfft::init (%i, {%i, %i, %i}, %i, {%i, %i, %i}, %i, *, *, %.9f)\n", 
			m_dim, 
			m_N[0], m_N[1], m_N[2],
			m_M,
			m_n[0], m_n[1], m_n[2],
			m,
			m_epsilon);

	for (int i = 0; i < NTHREADS || i < m_Nc; i++)
		nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan[i], &m_iplan[i]);
	// --------------------------------------

	m_initialised = true;

	printf ("... done.\n\n");

	return error;

}


RRSModule::error_code
CGSENSE::Prepare () {

	RRSModule::error_code error = OK;

	Matrix<cxfl>&   sens    = GetCXFL("sens");
	Matrix<double>& weights = GetRLDB("weights");
	Matrix<double>& kspace  = GetRLDB("kspace");

	m_ncs = new NCSENSE<cxdb> (sens, m_M, 1.0e-6, 20);

	m_ncs->KSpace (GetRLDB ("kspace"));
	m_ncs->Weights (GetRLDB ("weights"));
	
	//FreeRLDB ("kspace");
	//FreeRLDB ("weights");
	//FreeCXFL ("sense")
	
	return error;



}


RRSModule::error_code
CGSENSE::Process () {

	RRSModule::error_code error = OK;

	Matrix<cxfl>&   data    = GetCXFL("data");
	Matrix<cxfl>&   sens    = GetCXFL("sens");
	Matrix<double>& weights = GetRLDB("weights");
	Matrix<double>& kspace  = GetRLDB("kspace");

	//image = ns ->* data;

	// CG matrices ----------------------------------------------------
	Matrix <cxfl>   p       = Matrix<cxfl>   (m_N[0],m_N[1],m_N[2]), q, r;

	// Intensity Correction -------------------------------------------
	m_intcor                = ones<double> (m_N[0],m_N[1],m_N[2]);

	// Set k-space and weights in FT plans and clear RAM --------------
	for (int i = 0; i < NTHREADS || i < m_Nc; i++) {

		memcpy (&(m_fplan[i].x[0]),  &kspace[0], m_fplan[i].d * m_fplan[i].M_total * sizeof(double));
		memcpy (&(m_iplan[i].w[0]), &weights[0],                m_fplan[i].M_total * sizeof(double));
		nfft::weights (&m_fplan[i], &m_iplan[i]);
		nfft::psi     (&m_fplan[i]);

	}

	FreeRLDB ("kspace");
	FreeRLDB ("weights");

	// Create test data if testcase (Incoming data is image space) ----
	if (m_testcase) {
		E  (phantom<cxfl>(m_N[0]), sens, m_intcor, m_fplan, data, m_dim);
		if (m_noise > 0.0)
			AddPseudoRandomNoise (data, (float) m_noise);
	} 

	IntensityCorrection (sens, m_intcor);
	
	// Out going images -----------------------------------------------
	Matrix<cxfl>& image = AddMatrix ("image", (Ptr<Matrix<cxfl> >) NEW (Matrix<cxfl>(m_N[0], m_N[1], m_N[2])));

	// Start CG routine and runtime -----------------------------------
	ticks cgstart = getticks();

	// First left side action -----------------------------------------
	EH (data, sens, m_intcor, m_fplan, m_iplan, m_epsilon, m_maxit, p, m_dim, false);
	r = p;
	q = p;

	int         iters = 0;

	if (m_verbose)
		memcpy (&image[iters*p.Size()], &p[0], p.Size() * sizeof(cxfl));

	// CG residuals storage and helper variables ----------------------
	std::vector<double> res;

	float       rn    = 0.0;
	float       rno   = 0.0;
	float       an    = 0.0;
	raw         rtmp  = raw(0.0,0.0);

    Matrix<cxfl> treg = eye<cxfl>(m_N[0]) * cxfl (m_lambda, 0);

	printf ("Processing CG-SENSE ...\n");

	// CG iterations (Pruessmann et al. (2001). MRM, 46(4), 638-51.) --

	an = pow(creal(Norm(p)), 2); 
	rn = an;

	Matrix<cxfl> a(m_N[0], m_N[1], m_N[2]);

	for (iters = 0; iters < m_cgmaxit; iters++) {

		res.push_back(rn/an);

		printf ("  %03i: CG residuum: %.9f\n", iters, res.at(iters));

		// Convergence ? ----------------------------------------------
		if (res.at(iters) <= m_cgeps) break;

		// EHE --------------------------------------------------------
		E  (p,    sens, m_intcor, m_fplan,                              data, m_dim);
		EH (data, sens, m_intcor, m_fplan, m_iplan, m_epsilon, m_maxit, q   , m_dim, false);
		q  += p * m_lambda;

		// Guess new gradient -----------------------------------------
		rtmp  = (rn / (p.dotc(q)));

		r    -= (q * rtmp);
		rno   = rn;
		rn    = pow(creal(Norm(r)), 2);
		if (std::isnan(rn))	break;
		
		a    += (p * rtmp);
		p    *= rn/rno;
		p    += r;

		// Verbose out put keeps all intermediate steps ---------------
		if (m_verbose) {

			if (iters)
				image.Expand(m_dim); 

			memcpy (&image[iters*a.Size()], &a[0], a.Size()*sizeof(cxfl));

		}
		
	}
	
	FreeCXFL ("sens");	
	
	// Report timimng -------------------------------------------------
	printf ("... done. WTime: %.4f seconds.\n\n", elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());

	// Verbose output needs to 
	if (m_verbose) {
		
		// CG residuals ------------------
		Matrix<double>& nrmse = AddMatrix ("nrmse", (Ptr<Matrix<double> >) NEW (Matrix<double> (iters,1)));
		memcpy (&nrmse[0], &res[0], iters * sizeof(double));
		
	} else {

		image = a;

		for (int i = 0; i < image.Size(); i++)
			image[i] *= m_intcor[i];

	}

	return error;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}




