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

#include "DirectMethod.hpp"

using namespace RRStrategy;

/**
 * @brief         Time reversal for RF 
 *
 * @param  signal Acquired signal
 * @param  jac    Jacobian determinant j(k) i.e. density compensation
 * @param  pulse  Excitation pulse(s)
 */
void TimeReverseRF (const Matrix<cplx>& signal, const Matrix<double>& jac, Matrix<cplx>& pulse) {
	
	size_t nt = signal.Dim(0);
	size_t nc = signal.Dim(1);
	
	for (size_t i = 0; i < nt; i++)
		for (size_t c = 0; c < nc; c++)
			pulse(i,c) = (signal(nt-1-i,c)*(float)jac[nt-1-i]);
	
}


/**
 * @brief         Time reversal for gradient trajectory
 * 
 * @param  acqgr  Acquisition gradients
 * @param  excgr  Excitation gradients
 */
void TimeReverseGR (const Matrix<double>& acqgr, Matrix<double>& excgr) {
	
	size_t nt = acqgr.Dim(1);
	
	for (size_t i = 0; i < nt; i++) {
		
		excgr(0,i) = -acqgr(0,nt-1-i); 
		excgr(1,i) = -acqgr(1,nt-1-i); 
		excgr(2,i) = -acqgr(2,nt-1-i);
		
	}
	
}


RRSModule::error_code 
DirectMethod::Init () {

    printf ("\nIntialising DirectMethod ...\n");
	
	m_verbose     = false;                  // verbose?
	m_np          = omp_get_num_threads();  // # procs
	m_dt          = 1e-6;                   // seconds 

	Attribute ("dt",      &m_dt);
    printf ("  delta t: %.6f \n", m_dt);

	Attribute ("verbose", (int*)&m_verbose);
    printf ("  verbose: %s \n", (m_verbose) ? "true" : "false");

	Attribute ("threads", &m_np);
    printf ("  # threads: %.i \n", m_np);

	Attribute ("ic", &m_ic);
    printf ("  intensity correction: %s \n", (m_ic) ? "true": "false");

	m_sim = new SimulationContext ();

	m_initialised = true;

    printf ("... done.\n");
	
	return RRSModule::OK;

}


RRSModule::error_code
DirectMethod::Finalise() {

	delete m_sim;
	
	FreeCplx("b1p");
	FreeCplx("b1m");
	FreeReal("ag");
	FreeReal("r");
	FreeReal("b0");
	FreeReal("target");
	FreeReal("sample");
	FreeReal("sb0");
	FreeReal("sr");
	FreeReal("j");

	return RRSModule::OK;

}


RRSModule::error_code
DirectMethod::Process     () { 

    printf ("Processing DirectMethod ...\n");
	ticks start = getticks();

	// Incoming
	Matrix<cplx>&   b1m    = GetCplx("b1m");
	Matrix<cplx>&   b1p    = GetCplx("b1p");

	Matrix<double>& ag     = GetReal("ag");
	Matrix<double>& r      = GetReal("r");
	Matrix<double>& b0     = GetReal("b0");
	Matrix<double>& target = GetReal("target");

	Matrix<double>& sample = GetReal("sample");
	Matrix<double>& sr     = GetReal("sr");
	Matrix<double>& sb0    = GetReal("sb0");
	Matrix<double>& j      = GetReal("j");
	
	// Outgoing
    Matrix<cplx>    res    = Matrix<cplx>  (ag.Dim(1),b1p.Dim(1));
    Matrix<cplx>&    rf    = AddCplx ("rf",      NEW (Matrix<cplx>  (ag.Dim(1),b1p.Dim(1))));
	Matrix<double>&   m    = AddReal ("magn",    NEW (Matrix<double>(        3, sr.Dim(1))));
	Matrix<double>   eg    = ag; 

	

	// Intensity correction (Vahedipour et al. MRM 2011)
	if (m_ic)
		IntensityCorrection (b1m, target);

	// Simulate Bloch receive mode
	m_sim->Simulate (b1p, b1m, rf, ag,  r, target,  b0, m_dt, ACQUIRE, m_verbose, m_np, res, m);

	// Time reversal
	TimeReverseRF (res, j, rf);
	TimeReverseGR (ag, eg);

	// Simulate Bloch transmit mode
	m_sim->Simulate (b1p, b1m, rf, eg, sr, sample, sb0, m_dt,  EXCITE, m_verbose, m_np, res, m);

	printf ("... done. Overall WTime: %.4f seconds.\n\n", elapsed(getticks(), start) / Toolbox::Instance()->ClockRate());

	return RRSModule::OK;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {

    return new DirectMethod;

}


extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {

	delete p;

}
