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

#include "KTPoints.hpp"
#include "PTXINIFile.hpp"
#include "Lapack.hpp"
#include "IO.hpp"

using namespace RRStrategy;


KTPoints::KTPoints  () {}


KTPoints::~KTPoints () {}


RRSModule::error_code
KTPoints::Init           ()  {

    printf ("Intialising KTPoints ...\n");
	
    RRSModule::error_code e = OK;
	
    // # kt points ---------------------------
    Attribute ("Nk",      &m_nk);
    printf ("  # kt-points: %i \n", m_nk);
	
    m_max_rf = (float*) malloc (100 * sizeof(float));
    for (int i = 0; i < m_nk; i++)
        m_max_rf [m_nk] = 0.0; 
	
    // # of spatial sites --------------------
    Attribute ("Ns",      &m_ns); 
    printf ("  # spatial sites: %i \n", m_ns);
	
    // # of transmit channels ----------------
    Attribute ("Nc",      &m_nc);
    printf ("  # transmitter: %i \n", m_nc);
    
    // rf pulse durations --------------------
    m_pd = (int*) malloc (100 * sizeof(int));
    Attribute ("pd",      &m_pd[0]);   
    printf ("  starting pulse durations: %ius \n", m_pd[0]*10);
	
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


RRSModule::error_code
KTPoints::Finalise       ()  {

    free (m_pd);
    free (m_max_rf);

    return OK;

}


RRSModule::error_code
KTPoints::Process        () { 

    printf ("Processing KTPoints ...\n");

    // On entry -------------------
    //
    // m_raw:     target pattern
    // m_kspace:  kt points
    // m_pixel:   b0 (Hz)
    // m_rhelper: b1+ maps
    // m_helper:  spatial positions
    // ----------------------------

    // On exit --------------------
    // 
    // m_raw:    Excitation profile
    // m_helper: RF and gradient pulses
    // ----------------------------

	Matrix<double>& k      = GetRLDB("k");
	Matrix<double>& r      = GetRLDB("r");
	Matrix<cxfl>&   b1     = GetCXFL("b1");
	Matrix<short>&  b0     = GetSHRT("b0");
	Matrix<cxfl>&   target = GetCXFL("target");

	printf ("k: %s\n", DimsToCString(k));

	m_ns = r.Dim(1);
	m_nk = k.Dim(1);
	m_nc = b1.Dim(1);

    for (int i = 1; i < m_nk; i++)
        m_pd[i] = m_pd[0];

    Matrix<cxfl>    solution;
    Matrix<cxfl>    tmp;
    Matrix<cxfl>    final;    
    Matrix<cxfl>    treg   =  Matrix<cxfl>::Id(m_nc * m_nk) * cxfl (m_lambda, 0);
	Matrix<cxfl>    ve;
	Matrix<cxfl>    vp;

	if (m_verbose) {
	    ve  = Matrix<cxfl>(m_ns,        m_maxiter);
		vp  = Matrix<cxfl>(m_nk * m_nc, m_maxiter);
	} else {
	    ve  = Matrix<cxfl>(m_ns,       1);
		vp  = Matrix<cxfl>(m_nk * m_nc,1);
	}

    bool        pulse_amp_ok = false;

    float       nrmse = 0.0;

	std::vector<double> res;
	int         gc    = 0;

    // Start clock ------------------------
    ticks vestart = getticks();
    
    while (!pulse_amp_ok) {
        
        Matrix<cxfl> m (m_ns, m_nk*m_nc);

        STA (k, r, b1, b0, m_nc, m_nk, m_ns, m_gd, m_pd, m);
		
        Matrix<cxfl> minv;

        minv  = m.prodt (m);
        minv += treg;
        minv  = Lapack::Pinv(minv);
        minv  = minv.prod (m, 'N', 'C');
        
        // Valriable exchange method --------------

        
        printf ("  Starting variable exchange method ...\n");

        for (int j = 0; j < m_maxiter; j++, gc++) {

            solution = minv->*(target);
            tmp      = m->*(solution);

            NRMSE (target, tmp, gc, nrmse);

            res.push_back (nrmse);
            
			if (m_verbose) memcpy (&ve(0,gc), &tmp.At(0), tmp.Size() * sizeof(cxfl)); 
			if (gc && m_breakearly && (res.at(gc) > res.at(gc-1) || res.at(gc) < m_conv)) break;
            
            final    = solution;

			PhaseCorrection (target, tmp);
            
        } 
        
        printf ("\n  ... done.\n  Checking pulse amplitudes: "); fflush(stdout);
        
        // Check max pulse amplitude -----------------

        RFLimits (final, m_pd, m_nk, m_nc, m_max_rf); 

        pulse_amp_ok = true;

        for (int i = 0; i < m_nk; i++)
            if (m_max_rf[i] > m_rflim) pulse_amp_ok = false;
        
        // Update Pulse durations if necessary -------

        if (!pulse_amp_ok) {
            
            printf ("Pulse amplitudes to high!\n  Updating pulse durations ... to "); fflush(stdout);
            
            for (int i = 0; i < m_nk; i++) {
                m_pd[i] = 1 + (int) (m_max_rf[i] * m_pd[i] / m_rflim); 
				printf ("%i ", 10*m_pd[i]); fflush(stdout);
			}
            
            printf ("[us] ... done.\n");
        } else 
			printf ("OK\n");
        
	} // End of pulse duration loop

    printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), vestart) / Toolbox::Instance()->ClockRate());
	
	// Put actual maximum RF amplitude into first cell

	for (int i = 1; i < m_nk; i++)
		if (m_max_rf[i] > m_max_rf[0])
			m_max_rf[0] = m_max_rf[i];
	// -----------------------------------

	// Assemble gradient and RF timing ---

	PTXTiming              (final, k, m_pd, m_gd, m_nk, m_nc, b1);
	// -----------------------------------

	// Write pulse file for Siemens sequences 
	// Assuming (Sagittal/Transversal A>>P) 

	stringstream ofname;
	
	ofname << m_ptxfname << ".sag_ap";
	PTXWriteSiemensINIFile (b1, 2, 3, m_nc, 10, m_max_rf[0], ofname.str(), "s");
	ofname.str("");
	ofname << m_ptxfname << ".tra_ap";
	PTXWriteSiemensINIFile (b1, 2, 3, m_nc, 10, m_max_rf[0], ofname.str(), "t");
	// -----------------------------------

	// Return NRMSE down the road --------

	r.Dim(COL) = gc;
	r.Dim(LIN) = 1;
	r.Reset();
	for (int i = 0; i < gc; i++)
		r[i] = res[i];
	// -----------------------------------

	// Excitation profile ----------------
	if (m_verbose) {
			
		target.Dim(1) = gc; 
	    target.Reset();
		memcpy (&target[0], &ve[0], gc * target.Dim(0) * sizeof(cxfl));

	} else
		
		target = tmp;
	// ----------------------------------

    return RRSModule::OK;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {

    return new KTPoints;

}



extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {

    delete p;

}
