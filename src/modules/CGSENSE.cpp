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
#include "SEM.hpp"
#include <math.h>

#include <vector>

std::string sides[3] = {"Nx", "Ny", "Nz"};

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
		nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan[i], &m_iplan[i], m_epsilon);
	// --------------------------------------

	m_initialised = true;

	printf ("... done.\n\n");

	return error;

}


RRSModule::error_code
CGSENSE::Process () {

	RRSModule::error_code error = OK;

	Matrix<cplx>&   data    = GetCplx("data");
	Matrix<cplx>&   sens    = GetCplx("sens");
	Matrix<double>& weights = GetReal("weights");
	Matrix<double>& kspace  = GetReal("kspace");
	
	// CG matrices ----------------------------------------------------
	Matrix <cplx> p = Matrix<cplx> (m_N[0],m_N[1],m_N[2]), q, r;

	// Add white noise? (Only for testing) ----------------------------
	if (m_noise > 0.0)
		AddPseudoRandomNoise (data, (float)m_noise);
	
	// Set k-space and weights in FT plans and clear RAM --------------
	for (int i = 0; i < NTHREADS || i < m_Nc; i++) {

		memcpy (&(m_fplan[i].x[0]),  &kspace[0], m_fplan[i].d * m_fplan[i].M_total * sizeof(double));
		memcpy (&(m_iplan[i].w[0]), &weights[0],                m_fplan[i].M_total * sizeof(double));
		nfft::weights (&m_fplan[i], &m_iplan[i]);
		nfft::psi     (&m_fplan[i]);

	}

	FreeReal ("kspace");
	FreeReal ("weights");

	// Create test data if testcase (Incoming data is image space) ----
	Matrix<cplx> stmp (m_M, m_Nc);
	if (m_testcase) 
		E  (data, sens, m_fplan, stmp, m_dim);
	else 
		stmp = data;

	FreeCplx("data");
	
	// Out going images -----------------------------------------------
	Matrix<cplx>&   image = AddCplx (  "image", NEW (Matrix<cplx>(m_N[0], m_N[1], m_N[2])));

	// Out going signals ----------------------------------------------
	Matrix<cplx>& signals = AddCplx ("signals", NEW (Matrix<cplx>(m_M, m_Nc)));

	// Start CG routine and runtime -----------------------------------
	ticks cgstart = getticks();

	// First left side action -----------------------------------------
	EH (stmp, sens, m_fplan, m_iplan, m_epsilon, m_maxit, p, m_dim, false);
	r = p;
	q = p;

	int         iters = 0;

	if (m_verbose) {
		memcpy (  &image[iters*   p.Size()],    &p[0],    p.Size() * sizeof(cplx));
		memcpy (&signals[iters*stmp.Size()], &stmp[0], stmp.Size() * sizeof(cplx));
	}

	// CG residuals storage and helper variables ----------------------
	std::vector<double> res;

	float       rn    = 0.0;
	float       an    = 0.0;
	raw         rtmp  = raw(0.0,0.0);

    Matrix<cplx> treg = Matrix<cplx>::Id(m_N[0]) * cplx (m_lambda, 0);

	printf ("Processing CG-SENSE ...\n");

	// CG iterations (Pruessmann et al. (2001). MRM, 46(4), 638-51.) --

	an = pow(p.Norm().real(), 2);
	Matrix<cplx> a(m_N[0], m_N[1], m_N[2]);
	Matrix<cplx> s(m_M, m_Nc);


	for (iters = 0; iters < m_cgmaxit; iters++) {

		rn = pow(r.Norm().real(), 2);

		res.push_back(rn/an);

		printf ("  %03i: CG residuum: %.9f\n", iters, res.at(iters));

		// Convergence ? ----------------------------------------------
		if (std::isnan(res.at(iters)) || res.at(iters) <= m_cgeps)
			break;

		// CG step ----------------------------------------------------
		E  (p,    sens, m_fplan,                              stmp, m_dim);
		EH (stmp, sens, m_fplan, m_iplan, m_epsilon, m_maxit, q   , m_dim, false);

		rtmp  = (rn / (p.dotc(q)));
		a    += (p    * rtmp);
		s    += (stmp * rtmp);
		r    -= (q    * rtmp);
		p    *= cplx(pow(r.Norm().real(), 2)/rn);
		p    += r;

		// Verbose out put keeps all intermediate steps ---------------
		if (m_verbose) {
			image.Expand(m_dim); 
			signals.Expand(2);
			memcpy (  &image[(iters+1)*a.Size()], &a[0], a.Size()*sizeof(cplx));
			memcpy (&signals[(iters+1)*s.Size()], &s[0], s.Size()*sizeof(cplx));
		}
		
	}
	
	FreeCplx ("sens");	
	
	// Report timimng -------------------------------------------------
	printf ("... done. WTime: %.4f seconds.\n\n", elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());

	// Verbose output needs to 
	if (m_verbose) {
		
		// CG residuals ------------------
		Matrix<double>& nrmse = AddReal ("nrmse", NEW (Matrix<double> (iters,1)));
		memcpy (&nrmse[0], &res[0], iters * sizeof(double));
		
	} else {

		image   = a;
		signals = s;

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




