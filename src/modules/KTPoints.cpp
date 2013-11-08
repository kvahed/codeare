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
#include "IOContext.hpp"

using namespace RRStrategy;

/**
 * @brief           Normalised root-means-squared error
 *
 * @param  target   Target magnetisation
 * @param  result   Achieved result
 * @param  iter     Iteration
 * @param  nrmse    Returned NRMSE
 */
inline float
NRMSE                         (Matrix<cxfl>& target, const Matrix<cxfl>& result) {

    float nrmse = 0.0;
    size_t i;

#pragma omp parallel for reduction(+:nrmse)
    for (size_t i = 0; i < numel(target); i++ )
        nrmse += pow(abs(target[i]) - abs(result[i]), 2);

    nrmse = sqrt(nrmse)/norm(target);
    
    return nrmse;

}


/**
 * @brief           Normalised root-means-squared error
 *
 * @param  target   Target magnetisation
 * @param  result   Achieved result
 * @param  iter     Iteration
 * @param  nrmse    Returned NRMSE
 */
inline float
NOHO                         (Matrix<cxfl>& target, const Matrix<cxfl>& result) {

    float nrmse = 0.0;
    size_t i;

#pragma omp parallel for default (shared) private (i) reduction(+:nrmse)
	for (i = 0; i < numel(target); i++ )
		nrmse += pow(abs(target[i]) - abs(result[i]), 2);
    
    return nrmse;

}



/**
 * @brief           Phase correction from off-resonance
 *
 * @param  target   Target magnetisation
 * @param  result   Achieved result
 * @return          Phase corrected 
 */
inline static void
PhaseCorrection (Matrix<cxfl>& target, const Matrix<cxfl>& result) {
    
    size_t n = target.Size();
    
#pragma omp parallel for
    for (size_t i = 0; i < n; i++) 
        target[i] = abs(target[i]) * result[i] / abs(result[i]);

}


/**
 * @brief           STA integral
 *
 * @param  ks       k
 * @param  r        r
 * @param  b1       b1 map
 * @param  b0       b0 map
 * @param  nc       # of transmit channels
 * @param  nk       # of kspace positions
 * @param  ns       # of spatial positions
 * @param  gd       Gradient duration
 * @param  pd       Pulse durations
 *
 * @return          STA matrix
 */
inline static Matrix<cxfl>
STA (const Matrix<float>& ks, const Matrix<float>& r, const Matrix<cxfl>& b1, const Matrix<float>& b0, 
     const size_t nc, const size_t nk, const size_t ns, const size_t gd, const Matrix<short>& pd, const size_t n) {

    container<float> d (nk);
	Matrix<cxfl>  m (ns,nc*nk);

    printf ("  Computing STA encoding matrix ..."); 
    fflush (stdout);
    
    for (size_t i = 0; i< nk; i++)
        d[i] = (i==0) ? pd[i] + gd : d[i-1] + pd[i] + gd;

    std::reverse (d.begin(), d.end());

    for (size_t i = 0; i < nk-1; i++)
        d[i] = 1.0e-5 * d[i+1] + 1.0e-5 * pd[i]/2;

    d[nk-1] = 1.0e-5 * pd[nk-1] / 2;

    cxfl pgd = cxfl(0, TWOPI * GAMMA * 10.0);

    #pragma omp for
            for (size_t k = 0; k < nk; k++) 
                for (size_t s = 0; s < ns; s++) {
                	cxfl eikr = pgd * exp (cxfl(0, ks(0,k)*r(0,s) + ks(1,k)*r(1,s) + ks(2,k)*r(2,s) + TWOPI * d[k] * b0(s)));
                	for (size_t c = 0; c < nc; c++)
                		m(s,c*nk+k) =  b1(s,c) * eikr;
                }
    
	printf ("... done. ");

	return m;

}


/**
 * @brief Construct actual pulses
 *
 *
 */
static inline void
PTXTiming (const Matrix<cxfl>& solution, const Matrix<float>& ks, const Matrix<short>& pd, const size_t gd, 
           const size_t nk, const size_t nc, Matrix<cxfl>& rf, Matrix<float>& grad) {
    
    
    // Total pulse duration --------------
    
    size_t tpd = 2;  // Start and end
    for (size_t i = 0; i < nk-1; i++) 
        tpd += (size_t) (pd[i] + gd);
    tpd += pd[nk-1];
    // -----------------------------------
    
    // Outgoing brepository ---------------
    
    rf   = Matrix<cxfl>  (tpd, nc);
    grad = Matrix<float> (tpd, 3);
    // -----------------------------------
    
    // RF Timing -------------------------

    for (size_t rc = 0; rc < nc; rc++) {
        
        size_t i = 1;
        
        for (size_t k = 0; k < nk; k++) {
            
            // RF action
            for (short p = 0; p < pd[k]; p++, i++)
                rf (i,rc) = solution(k + rc*nk) / (float) pd[k] * cxfl(1000.0,0.0);

            // Gradient action, no RF
            if (k < nk-1)
                for (size_t g = 0; g <    gd; g++, i++)
                    grad (i,0) = 0.0;

        }
        
    }
    // -----------------------------------
    
    // Gradient and slew -----------------
    
    float gr = 0.0;
    float sr = 0.0;

    for (size_t gc = 0; gc < 3; gc++) {
        
        size_t i = 1;
        
        for (size_t k = 0; k < nk-1; k++) {
            
            // RF action, no gradients
            for (short p = 0; p < pd[k]; p++, i++)
                grad (i,gc) = 0.0; 
            
            if (k < nk-1) 
                
                for (size_t g = 0; g <    gd; g++, i++) {
                    
                    sr = (k+1 < nk) ? ks(gc,k+1) - ks(gc,k) : - ks(gc,i);
                    sr = 4.0 * sr / (2 * PI * GAMMA * 1e6 * gd * gd * 1.0e-5 * 1.0e-5);
                    sr =       sr / 100.0;
                    
                    // Gradient action 
                    if(g < gd/2)             // ramp up
                        gr = sr * (0.5 + g);
                    else if (g < gd/2+1)     // flat top
                        0;
                    else                     // ramp down
                        gr -= sr;

                    grad (i,gc) = gr; 
                }
            
        } 
        
    }
    // ----------------------------------
    
}



static inline void
KTPSolve (const Matrix<cxfl>& m, Matrix<cxfl>& target, Matrix<cxfl>& final,
          Matrix<cxfl>& solution, const double& lambda, const size_t& mxit, 
          const float& conv, const bool& breakearly, size_t& gc, 
		  Matrix<float>& res) {

    ticks start = getticks();
    printf ("  Starting variable exchange method ...\n");
    
    Matrix<cxfl> treg = lambda *  eye<cxfl>(size(m,1));
    Matrix<cxfl> minv;

    minv  = m.prodt(m); 
    minv += treg;
    minv  = inv(minv);
    minv  = minv.prod (m, 'N', 'C');
    
    // Valriable exchange method --------------
    for (size_t j = 0; gc < mxit; j++, gc++) {
		
        solution = minv ->* target;
        final    = m    ->* solution;
        
        res[gc]  = NRMSE (target, final);
        PhaseCorrection  (target, final);
		
		if (j % 5 == 0 && j > 0)
			printf ("\n");

		printf ("    %04zu %.6f", gc, res[gc]); 

		fflush (stdout);

        if ((gc && j && breakearly && res[gc] > res[gc-1]) || res[gc] < conv || isnan(res[gc])) {
			gc++;
            break; 
		} 
        
    } 
    
    printf ("\n  ... done. WTime: %.4f seconds.\n", elapsed(getticks(), start) / Toolbox::Instance()->ClockRate());
    
}



inline static bool 
CheckAmps (const Matrix<cxfl>& solution, Matrix<short>& pd, const size_t& nk, 
		   const size_t& nc, Matrix<float>& max_rf, const float& rflim) {

	bool amps_ok = true;

	printf ("  Checking pulse amplitudes: "); 
	fflush(stdout);
	
    for (size_t i = 0; i < nk; i++) {

        max_rf[i] = 0.0;
        
        for (size_t j = 0; j < nc; j++)
            if (max_rf[i] < abs (solution[i+nk*j]) / (float)(10.0*pd[i])) 
                max_rf[i] = abs (solution[i+nk*j]) / (float)(10.0*pd[i]);

        max_rf[i] *= 100.0;

    }
    
	for (size_t i = 0; i < nk; i++)
		if (max_rf[i] > rflim) amps_ok = false;
	
	// Update Pulse durations if necessary -------
    
	if (!amps_ok) {
		
		printf ("Pulse amplitudes to high!\n  Updating pulse durations ... to "); fflush(stdout);
        
		for (size_t i = 0; i < nk; i++) {
			pd[i] = 1 + (short) (max_rf[i] * pd[i] / rflim); 
			printf ("%i ", 10*pd[i]); fflush(stdout);
		}
        
		printf ("[us] ... done.\n\n");
		
	} else 

		printf ("OK\n\n");

	return amps_ok;
	
}


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
	pd *= 5;

    Matrix<cxfl>    solution; // Solution for timing calculation
    Matrix<cxfl>&   final = AddMatrix<cxfl> ("final");   // Excitation profile

    Matrix<cxfl>&   rf    = AddMatrix<cxfl> ("rf"); 
    Matrix<float>&  grad  = AddMatrix<float> ("grad");

    Matrix<cxfl> m;// (ns, nk*nc); // STA system encoding matrix

    bool        amps_ok = false;
    size_t      gc    = 0;      // Global counter for VE iterations

	Matrix<float>& res  = AddMatrix<float> ("nrmse");
    res = Matrix<float> (m_maxiter,1);

    SimpleTimer st ("KT-Points");

    while (!amps_ok) {
        
		// Compute SEM
        m = STA (k, r, b1, b0, nc, nk, ns, m_gd, pd, gc);

		// Solve KTPoints
        KTPSolve (m, target, final, solution, m_lambda, m_maxiter, m_conv, m_breakearly, gc, res);
		if (isnan(res(gc)))
			break;
    
        // Check max pulse amplitude -----------------        
		amps_ok = CheckAmps(solution, pd, nk, nc, max_rf, m_rflim);
        
    } // End of pulse duration loop

    // Put actual maximum RF amplitude into first cell
    for (size_t i = 1; i < nk; i++)
        if (max_rf[i] > max_rf[0])
            max_rf[0] = max_rf[i];
    // -----------------------------------
    // Assemble gradient and RF timing ---

    PTXTiming (solution, k, pd, m_gd, nk, nc, rf, grad);
    // -----------------------------------

    // Write pulse file for Siemens sequences 
    // Assuming (Sagittal/Transversal A>>P) 

    stringstream ofname;

	// Resize NRMSE vector
	res = resize (res,gc,1);

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
