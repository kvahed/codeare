/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
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

#include "SpatialDomain.hpp"
#include "PTXINIFile.hpp"

using namespace RRStrategy;


SpatialDomain::SpatialDomain  () {}


SpatialDomain::~SpatialDomain () {}


RRSModule::error_code
SpatialDomain::Init           ()  {

    printf ("Intialising SpatialDomain ...\n");
	
    RRSModule::error_code e = OK;
	
    // # kt points ---------------------------
    Attribute ("Nk",      &m_nk);
    printf ("  # kt-points: %i \n", m_nk);
	
    m_max_rf = (float*) malloc (1000 * sizeof(float));
    for (int i = 0; i < m_nk; i++)
        m_max_rf [m_nk] = 0.0; 
	
    // # of spatial sites --------------------
    Attribute ("Ns",      &m_ns); 
    printf ("  # spatial sites: %i \n", m_ns);
	
    // # of transmit channels ----------------
    Attribute ("Nc",      &m_nc);
    printf ("  # transmitter: %i \n", m_nc);
    
    // rf pulse durations --------------------
    m_pd = (int*) malloc (m_nk * sizeof(int));
    Attribute ("pd",      &m_pd[0]);   
    for (int i = 1; i < m_nk; i++)
        m_pd[i] = m_pd[0];
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

	// ----------------------------------------
	
    printf ("... done.\n\n");
	
    return e;

}


RRSModule::error_code
SpatialDomain::Finalise       ()  {

    free (m_pd);
    free (m_max_rf);

    return OK;

}


RRSModule::error_code
SpatialDomain::Process        () { 

    printf ("Processing SpatialDomain ...\n");

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

    Matrix<raw> solution;
    Matrix<raw> tmp;
    Matrix<raw> final;    
    Matrix<raw> treg         =  Matrix<raw>::Id(m_nc * m_nk) * raw (m_lambda, 0);
	Matrix<raw> ve;
	Matrix<raw> vp;

	if (m_verbose) {
	    ve  = Matrix<raw>(m_ns,      m_maxiter);
		vp  = Matrix<raw>(m_nk*m_nc, m_maxiter);
	} else {
	    ve  = Matrix<raw>(m_ns,     1);
		vp  = Matrix<raw>(m_nk*m_nc,1);
	}

    bool        pulse_amp_ok = false;

    float       nrmse = 0.0;

	std::vector<double> res;
	int         gc    = 0;

    // Start clock ------------------------
    ticks vestart = getticks();
    
    while (!pulse_amp_ok) {
        
        Matrix<raw> m (m_ns, m_nk*m_nc);

        STA (&m_kspace, &m_helper, &m_rhelper, &m_pixel, m_nc, m_nk, m_ns, m_gd, m_pd, &m);
        Matrix<raw> minv = m.tr();

        minv  = minv->*(m);
        minv  = minv + treg;
        minv  = minv.Pinv();
        minv  = minv.prodt(m);
        
        // Valriable exchange method --------------

        
        printf ("  Starting variable exchange method ...\n");

        for (int j = 0; j < m_maxiter; j++, gc++) {

            solution = minv->*(m_raw);
            tmp      = m->*(solution);

            NRMSE (&m_raw, &tmp, gc, &nrmse);
            res.push_back (nrmse);
            
			if (m_verbose)
				memcpy (&ve.At(0,gc), &tmp.At(0), tmp.Size() * sizeof(raw)); 

            if (gc > 0 && (res.at(gc) > res.at(gc-1) || res.at(gc) < m_conv)) 
				break;
            
            final    = solution;
            PhaseCorrection (&m_raw, &tmp);
            
        } 
        
        printf ("\n... done. Checking pulse amplitudes ... \n");
        
        // Check max pulse amplitude -----------------
        RFLimits (&final, m_pd, m_nk, m_nc, m_max_rf); 

        pulse_amp_ok = true;
        
        for (int i = 0; i < m_nk; i++) 
            if (m_max_rf[i] > m_rflim)
                pulse_amp_ok = false;
        
        // Update Pulse durations if necessary -------
        if (!pulse_amp_ok) {
            
            printf ("... done. Pulse amplitudes to high! Updating pulse durations ...\n");
            
            for(int i=0; i < m_nk;i++)
                m_pd[i] = 1 + (int) (m_max_rf[i] * m_pd[i] / m_rflim);
            
            printf ("... done\n.");
            
        } 
        
    } // End of pulse duration loop

    printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), vestart) / ClockRate());
	
	// Put actual maximum RF amplitude into first cell

	for (int i = 1; i < m_nk; i++)
		if (m_max_rf[i] > m_max_rf[0])
			m_max_rf[0] = m_max_rf[i];
	// -----------------------------------

	// Assemble gradient and RF timing ---

	PTXTiming              (&final, &m_kspace, m_pd, m_gd, m_nk, m_nc, &m_rhelper);
	// -----------------------------------
	

	// Write pulse file for Siemens sequences

	PTXWriteSiemensINIFile (&m_rhelper, 3, 3, m_nc, 10, m_max_rf[0], &m_ptxfname);
	// -----------------------------------

	// Return NRMSE down the road --------

	m_helper.Dim(COL) = gc;
	m_helper.Dim(LIN) = 1;
	m_helper.Reset();
	for (int i = 0; i < gc; i++)
		m_helper.At(i) = res.at(i);
	// -----------------------------------

	// Excitation profile ----------------
	if (m_verbose) {
			
		m_raw.Dim(1) = gc; 
		m_raw.Reset();
		memcpy (&m_raw.At(0), &ve.At(0), gc * m_raw.Dim(0) * sizeof(raw));

	} else
		
		m_raw = tmp;
	// ----------------------------------
	
    return RRSModule::OK;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {

    return new SpatialDomain;

}



extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {

    delete p;

}
