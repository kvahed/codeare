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

    Attribute ("dt",      &m_dt);
    printf ("  delta t: %.6f \n", m_dt);

    Attribute ("verbose", (int*)&m_verbose);
    printf ("  verbose: %s \n", (m_verbose) ? "true" : "false");

    Attribute ("threads", &m_np);
    printf ("  # threads: %.i \n", m_np);

    Attribute ("ic", &m_ic);
    printf ("  intensity correction: %s \n", (m_ic) ? "true": "false");

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
    if (m_ic) IntensityCorrection (GetCplx("b1m"), GetFloat("target"));

    sb.tb1  = m_cplx["b1m"];
    sb.sb1  = m_cplx["b1p"];

    sb.agr  = m_float["ag"];
    sb.tm   = m_float["target"];
    sb.sm   = m_float["sample"];
    sb.tb0  = m_float["b0"];
    sb.sb0  = m_float["sb0"];
    sb.tr   = m_float["r"];
    sb.sr   = m_float["sr"];
    sb.jac  = m_float["j"];

    sb.np   = m_np;
    sb.mode = m_mode;
    sb.dt   = m_dt;
    sb.v    = m_verbose;

    // Outgoing
    AddCplx  (  "rf", sb.rf   = NEW (Matrix<cplx>  (sb.agr->Dim(1), sb.tb1->Dim(1))));
    AddFloat ("magn", sb.magn = NEW (Matrix<float> (             3,  sb.sr->Dim(1))));

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
