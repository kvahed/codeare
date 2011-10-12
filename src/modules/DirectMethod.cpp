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
#include "BlochSim.hpp"

using namespace RRStrategy;

RRSModule::error_code 
DirectMethod::Init () {

	m_verbose     = false;                  // verbose?
	m_np          = omp_get_num_threads();  // # procs
	m_dt          = 1e-6;                   // seconds 

	Attribute ("dt",      &m_dt);
	Attribute ("verbose", (int*)&m_verbose);
	Attribute ("threads", &m_np);

	m_initialised = true;

	return RRSModule::OK;

}



RRSModule::error_code
DirectMethod::Process     () { 

	// Nomen est omen
	Matrix<double>& k      = GetReal("k");
	Matrix<double>& r      = GetReal("r");
	Matrix<double>& b0     = GetReal("b0");
	Matrix<cplx>&   b1     = GetCplx("b1");
	Matrix<double>& target = GetReal("target");

	// Dummy: We want to use the complex conjugated transmit maps and receive maps 
	Matrix<cplx>    b1p;
	Matrix<cplx>    rf;
	Matrix<double>  m;

	// Complex conjugate transmit sensitivities for receive.
	b1.Conj();

	// Resulting signal
    Matrix<cplx>&   res    = AddCplx ("rxm",  NEW (Matrix<cplx>()));

	// Simulate Bloch
	Simulate (b1p, b1, rf, k, r, target, b0, m_dt, false, m_verbose, m_np, res, m);

	return RRSModule::OK;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {

    return new DirectMethod;

}


extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {

	delete p;

}
