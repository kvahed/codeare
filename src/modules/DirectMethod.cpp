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

#include "DirectMethod.hpp"
#include "Interpolate.hpp"

using namespace RRStrategy;

error_code 
DirectMethod::Init () {

    printf ("\nIntialising DirectMethod ...\n");
    
    m_verbose     = false;                  // verbose?
    m_np          = omp_get_num_threads();  // # procs
    m_dt          = 1e-6;                   // seconds 
	m_np          = 8;

    Attribute ("dt",      &m_dt);
    printf ("  delta t: %.6f \n", m_dt);

    Attribute ("cgeps",      &m_cgeps);
    printf ("  CG eps: %.6f \n", m_cgeps);

    Attribute ("cgiter",      &m_cgiter);
    printf ("  CG iters: %.i \n", m_cgiter);

    Attribute ("lambda",      &m_lambda);
    printf ("  Tikhonov factor: %f \n", m_lambda);

    Attribute ("verbose", &m_verbose);
    printf ("  verbose: %s \n", (m_verbose) ? "true" : "false");

    Attribute ("cb0", &m_cb0);
    printf ("  cb0: %s \n", (m_cb0) ? "true" : "false");

    Attribute ("threads", &m_np);
    printf ("  # threads: %.i \n", m_np);

    Attribute ("mode", &m_mode);
    printf ("  mode: %s \n", (m_mode) ? "true": "false");

    m_initialised = true;

    printf ("... done.\n");


    return OK;

}


error_code
DirectMethod::Finalise() {

	ReconStrategy::Finalise();
    return OK;

}


error_code
DirectMethod::Process     () { 

    printf ("Processing DirectMethod ..."); fflush (stdout);

    SimulationBundle sb;
	
    ticks           start  = getticks();

    
    Matrix<float>& t = Get<float> ("t");
	
	float t0, tT;
	t0 = 0.0;
	tT = t [numel(t)-1];
	Matrix<float> nt = linspace(t0, tT, tT/m_dt);
	
	Matrix<float> ngx = !Get<float> ("g");
	Matrix<float> nj  =  Get<float> ("j");
	
	ngx = !interp1 (t, ngx, nt, INTERP::LINEAR);
	nj  = interp1  (t, nj, nt);
	sb.g     = boost::make_shared<Matrix<float> >(ngx);

    sb.b1    = boost::make_shared<Matrix<cxfl> >(Get<cxfl> ("b1"));
    sb.tmxy  = boost::make_shared<Matrix<cxfl> >(Get<cxfl> ("tmxy"));
    sb.smxy  = boost::make_shared<Matrix<cxfl> >(Get<cxfl> ("smxy"));

    sb.tmz   = boost::make_shared<Matrix<float> >(Get<float> ("tmz"));
    sb.smz   = boost::make_shared<Matrix<float> >(Get<float> ("smz"));
    sb.roi   = boost::make_shared<Matrix<float> >(Get<float> ("roi"));
    sb.b0    = boost::make_shared<Matrix<float> >(Get<float> ("b0"));
    sb.r     = boost::make_shared<Matrix<float> >(Get<float> ("r"));
    sb.jac   = boost::make_shared<Matrix<float> >(nj);

    sb.np     = m_np;
    sb.mode   = m_mode;
    sb.dt     = m_dt;
    sb.v      = m_verbose;
	sb.cgeps  = m_cgeps;
	sb.cgit   = m_cgiter;
	sb.lambda = m_lambda;
	sb.cb0    = m_cb0;

    // Outgoing

    Matrix<cxfl> rf (sb.g->Dim(1), sb.b1->Dim(1));
    sb.rf = boost::make_shared<Matrix<cxfl> >(rf);
    Matrix<cxfl> mxy (sb.r->Dim(1), 1);
    sb.mxy = boost::make_shared<Matrix<cxfl> >(mxy);
    Matrix<float> mz (sb.r->Dim(1), 1);
    sb.mz = boost::make_shared<Matrix<float> >(mz);

    // Initialise CPU/GPU simulator
    SimulationContext sc (&sb);

    // Simulate
	printf (" simulating ...\n"); fflush(stdout);
    sc.Simulate();

    printf ("... done. Overall WTime: %.4f seconds.\n\n", elapsed(getticks(), start) / Toolbox::Instance()->ClockRate());
    return OK;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {

    return new DirectMethod;

}


extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {

    delete p;

}
