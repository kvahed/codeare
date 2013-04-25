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

#include "KTPoints.hpp"
#include "PTXINIFile.hpp"

using namespace RRStrategy;


KTPoints::KTPoints  () :
        m_verbose   (false),
        m_rflim     (1.0),
        m_conv      (1.0e-6),
        m_lambda    (1.0e-6),
        m_breakearly(true),
        m_gd        (1.0e-5),
        m_max_rf    (0),
        m_maxiter   (1000),
        ns(0), nk(0), nc(0) {

}


KTPoints::~KTPoints ()     {}


error_code
KTPoints::Init      ()     {

    printf ("Intialising KTPoints ...\n");
    
    error_code e = OK;
    
    // rf pulse durations --------------------
    /*
    m_pd = (int*) malloc (sizeof(int));
    Attribute ("pd",      &m_pd[0]);   
    printf ("  starting pulse durations: %ius \n", m_pd[0]*10);
    */
    // gradient pulse duration ----------------
    Attribute ("gd",      &m_gd);      
    printf ("  gradient blip durations: %ius \n", m_gd*10);
    
    // Max # of iterations --------------------
    Attribute ("maxiter", &m_maxiter);
    printf ("  maximum iterations: %i \n", m_maxiter);
    
    // # Tikhonov parameter from L-curve learning 
    // ca. [0.005...0.05] for human -----------
    Attribute ("lambda",  &m_lambda);  
    printf ("  tikhonov parameter: %.4f \n", m_lambda);
    
    // # RF limit (usualy in [0 1]) -----------
    Attribute ("rflim",   &m_rflim);  
    printf ("  RF amplitude limit: %.4f \n", m_rflim);
    
    // Convergence limit ----------------------
    Attribute ("conv",    &m_conv); 
    printf ("  Convergence criterium: %.4f \n", m_conv);
    
    // Pulse orientation ----------------------
    m_orient = Attribute ("orientation");
    printf ("  orientation: %s \n", m_orient.c_str());

    // PTX file name --------------------------
    m_ptxfname = Attribute ("ptxfname");
    printf ("  ptx file name: %s \n", m_ptxfname.c_str());

    // PTX file name --------------------------
    Attribute ("verbose", &m_verbose);
    printf ("  verbosity: %i \n", m_verbose);

    // Break loop early --------------------------
    Attribute ("breakearly", &m_breakearly);
    printf ("  break early: %i \n", m_breakearly);

    // ----------------------------------------
    
    printf ("... done.\n\n");
    
    return e;

}


error_code
KTPoints::Finalise  ()     {

    return OK;

}


error_code
KTPoints::Process   ()     {

    printf ("Processing KTPoints ...\n");

    // On entry -------------------
    //
    // target: target pattern
    // k:      kt points
    // b0:     b0 (Hz)
    // b1:     b1+ maps
    // r:      spatial positions
    // ----------------------------

    // On exit --------------------
    // 
    // rf:     RF pulses
    // grad:   Gradient pulses
    // ----------------------------

    Matrix<float>& k      = Get<float> ("k");
    Matrix<float>& r      = Get<float> ("r");
    Matrix<cxfl>&  b1     = Get<cxfl>  ("b1");
    Matrix<float>& b0     = Get<float> ("b0");
    Matrix<cxfl>&  target = Get<cxfl>  ("target");

    size_t ns = size( r,1); // # of spatial positions
    size_t nk = size( k,1); // # of kt points
    size_t nc = size(b1,1); // # of RF channels

    printf ("  # spatial sites: %zu \n", ns);
    printf ("  # transmitter: %zu \n", nc);
    printf ("  # kt points: %zu \n", nk);

    Matrix<float> max_rf (nk,1);
    Matrix<short> pd = ones<short>(nk,1); // Starting with shortest pulses possible 

    Matrix<cxfl>    solution; // Solution for timing calculation
    Matrix<cxfl>    final;    // Excitation profile

    Matrix<cxfl>&   rf     = AddMatrix (  "rf",  (Ptr<Matrix<cxfl> >)  NEW (Matrix<cxfl>  ())); 
    Matrix<float>&  grad   = AddMatrix ("grad",  (Ptr<Matrix<float> >) NEW (Matrix<float> ()));

    Matrix<cxfl> m (ns, nk*nc); // STA system encoding matrix

    bool        amps_ok = false;
    size_t      gc    = 0;      // Global counter for VE iterations

	Matrix<float>& res     = AddMatrix ("nrmse",  (Ptr<Matrix<float> >) NEW (Matrix<float> (m_maxiter,1)));

    ticks vestart = getticks();
    printf ("Starting KT-Points algorithm ...\n");
    
    while (!amps_ok) {
        
		// Compute SEM
        STA   (k, r, b1, b0, nc, nk, ns, m_gd, pd, m);

		// Solve KTPoints
        KTPSolve (m, target, final, solution, m_lambda, m_maxiter, m_conv, m_breakearly, gc, res);
    
        // Check max pulse amplitude -----------------        
		amps_ok = CheckAmps(solution, pd, nk, nc, max_rf, m_rflim);
        
    } // End of pulse duration loop

    printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), vestart) /
            Toolbox::Instance()->ClockRate());
    
    // Put actual maximum RF amplitude into first cell
    for (int i = 1; i < nk; i++)
        if (max_rf[i] > max_rf[0])
            max_rf[0] = max_rf[i];
    // -----------------------------------
    // Assemble gradient and RF timing ---

    PTXTiming (solution, k, pd, m_gd, nk, nc, rf, grad);
    // -----------------------------------

    // Write pulse file for Siemens sequences 
    // Assuming (Sagittal/Transversal A>>P) 

    stringstream ofname;

    ofname << m_ptxfname << ".sag_ap";
    PTXWriteSiemensINIFile (rf, grad, 2, 3, nc, 10, max_rf[0], ofname.str(), "s");
    ofname.str("");
    ofname << m_ptxfname << ".tra_ap";
    PTXWriteSiemensINIFile (rf, grad, 2, 3, nc, 10, max_rf[0], ofname.str(), "t");
    // -----------------------------------
    return OK;

}


// Factory
extern "C" DLLEXPORT ReconStrategy*
create  ()                 {

    return new KTPoints;

}


// Trash bin
extern "C" DLLEXPORT void
destroy (ReconStrategy* p) {

    delete p;

}
