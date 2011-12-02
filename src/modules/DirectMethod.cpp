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

RRSModule::error_code 
DirectMethod::Init () {

    printf ("\nIntialising DirectMethod ...\n");
    
    m_verbose     = false;                  // verbose?
    m_np          = omp_get_num_threads();  // # procs
    m_dt          = 1e-6;                   // seconds 
	m_np          = 8;

    Attribute ("dt",      &m_dt);
    printf ("  delta t: %.6f \n", m_dt);

    Attribute ("rfsc",      &m_rfsc);
    printf ("  rf scaling factor: %.6f \n", m_rfsc);

    Attribute ("cgeps",      &m_cgeps);
    printf ("  CG eps: %.6f \n", m_cgeps);

    Attribute ("cgiter",      &m_cgiter);
    printf ("  CG iters: %.i \n", m_cgiter);

    Attribute ("verbose", &m_verbose);
    printf ("  verbose: %s \n", (m_verbose) ? "true" : "false");

    Attribute ("threads", &m_np);
    printf ("  # threads: %.i \n", m_np);

    Attribute ("ic", &m_ic);
    printf ("  intensity correction: %s \n", (m_ic) ? "true": "false");

    Attribute ("mode", &m_mode);
    printf ("  mode: %s \n", (m_mode) ? "true": "false");

    m_initialised = true;

    printf ("... done.\n");
    
    return RRSModule::OK;

}


RRSModule::error_code
DirectMethod::Finalise() {

	ReconStrategy::Finalise();
    return RRSModule::OK;

}


RRSModule::error_code
DirectMethod::Process     () { 

    printf ("Processing DirectMethod ...\n");

    SimulationBundle sb;
    ticks           start  = getticks();

    // Intensity correction for single run
    if (m_ic) 
		IntensityCorrection (GetCplx("b1"), GetCplx("tmxy"), GetFloat("tmz"));

    sb.b1    = m_cplx ["b1"];
    sb.tmxy  = m_cplx ["tmxy"];
    sb.smxy  = m_cplx ["smxy"];

    sb.agr   = m_float["ag"];
    sb.tmz   = m_float["tmz"];
    sb.smz   = m_float["smz"];
    sb.roi   = m_float["roi"];
    sb.b0    = m_float["b0"];
    sb.r     = m_float["r"];
    sb.jac   = m_float["j"];

    sb.np    = m_np;
    sb.mode  = m_mode;
    sb.dt    = m_dt;
    sb.rfsc  = m_rfsc;
    sb.v     = m_verbose;
	sb.cgeps = m_cgeps;
	sb.cgit  = m_cgiter;

    // Outgoing
    AddCplx  ( "rf", sb.rf  = NEW (Matrix<cplx>  (sb.agr->Dim(1), sb.b1->Dim(1))));
    AddCplx  ("mxy", sb.mxy = NEW (Matrix<cplx>  (             1, sb.r->Dim(1))));
    AddFloat ( "mz", sb.mz  = NEW (Matrix<float> (             1, sb.r->Dim(1))));

    // Initialise CPU/GPU simulator
    SimulationContext sc (&sb);

    // Simulate
    sc.Simulate();

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
