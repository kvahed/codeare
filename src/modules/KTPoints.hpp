/*
 *  codeare Copyright (C) 2010-2012 Kaveh Vahedipour
 *                                  Forschungszentrum Juelich, Germany
 *  
 *  Original code stolen from Martijn Cloos
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

#ifndef __KT_POINTS_HPP__
#define __KT_POINTS_HPP__

#include "ReconStrategy.hpp"
#include "Toolbox.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "IO.hpp"
#include "linalg/Lapack.hpp"

/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

    /**
     * @brief Spatial domain method for RF pulse generation with variable exchange method<br/>
     *        Cloos MA, Boulant N, Luong M, Ferrand G, Giacomini E, Le Bihan D, Amadon A.<br/>
     *        kT-points: short three-dimensional tailored RF pulses for flip-angle homogenization over an extended volume. MRM:2012; 67(1), 72���80.
     */
    class KTPoints : public ReconStrategy {
        
        
    public:
        
        /**
         * @brief Default constructor
         */
        KTPoints  ();
        
        /**
         * @brief Default destructor
         */
        virtual 
        ~KTPoints ();
        
        
        /**
         * @brief Process
         */
        virtual error_code
        Process ();

        
        /**
         * @brief Initialise
         */
        virtual error_code
        Init ();
        

        /**
         * @brief Finalise
         */
        virtual error_code
        Finalise ();



    private:
 
        int           nc;       /**< @brief Transmit channels       */
        Matrix<short> m_pd;       /**< @brief Pulse durations         */
        int           m_gd;       /**< @brief Gradient durations      */
        int           ns;       /**< @brief # Spatial positions     */
        int           nk;       /**< @brief # kt-points             */
        int           m_maxiter;  /**< @brief # Variable exchange method iterations */
        int           m_verbose;  /**< @brief Verbose output. All intermediate results. */
        int           m_breakearly;  /**< @brief Break search with first diverging step */

        double        m_lambda;   /**< @brief Tikhonov parameter      */
        double        m_rflim;    /**< @brief Maximum rf amplitude    */
        double        m_conv;     /**< @brief Convergence criterium   */
        
        float*        m_max_rf;   /**< @brief Maximum reached RF amps */

        std::string   m_orient;   /**< @brief Orientation             */ 
        std::string   m_ptxfname; /**< @brief PTX file name           */

    };

}
#endif /* __DUMMY_RECON_H__ */


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

    for (size_t i = 0; i < numel(target); i++ )
        nrmse = nrmse + pow(abs(target[i]) - abs(result[i]), 2);
    
    nrmse = sqrt(nrmse)/norm(target);
    
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
 * @param  m        Out: m_xy
 */
inline static void
STA (const Matrix<float>& ks, const Matrix<float>& r, const Matrix<cxfl>& b1, const Matrix<float>& b0, 
     const int& nc, const int& nk, const int& ns, const int gd, const Matrix<short>& pd, Matrix<cxfl>& m) {

    ticks start = getticks();
    printf ("  Computing STA encoding matrix ..."); 
    fflush (stdout);
    
    vector<float> d (nk);
    vector<float> t (nk);

    for (int i = 0; i< nk; i++)
        d[i] = (i==0) ? pd[i] + gd : d[i-1] + pd[i] + gd;

    for (int i = 0; i< nk; i++)        
        t[i] = d [nk-i-1];
    
    for (int i = 0; i < nk-1; i++)
        d[i] = 1.0e-5 * t[i+1] + 1.0e-5 * pd[i]/2;

    d[nk-1] = 1.0e-5 * pd[nk-1] / 2;

    float pgd = 2.0 * PI * 4.2576e7 * 1.0e-5; 

#pragma omp parallel default (shared) 
    {

#pragma omp for
        for (int c = 0; c < nc; c++) 
            for (int k = 0; k < nk; k++) 
                for (int s = 0; s < ns; s++) 
                    m(s, c*nk + k) = 
                        // b1 (s,c)
                        pgd * b1(s,c) *
                        // off resonance: exp (2i\pidb0dt)  
                        exp (cxfl(0, 2.0 * PI * d[k] * (float) b0(s))) *
                        // encoding: exp (i k(t) r)
                        exp (cxfl(0, (ks(0,k)*r(0,s) + ks(1,k)*r(1,s) + ks(2,k)*r(2,s))));
        
    }

    printf (" done. WTime: %.4f seconds.\n", elapsed(getticks(), start) / Toolbox::Instance()->ClockRate());

}


/**
 * @brief Construct actual pulses
 *
 *
 */
static inline void
PTXTiming (const Matrix<cxfl>& solution, const Matrix<float>& ks, const Matrix<short>& pd, const int& gd, 
           const int& nk, const int& nc, Matrix<cxfl>& rf, Matrix<float>& grad) {
    

    // Total pulse duration --------------
    
    int tpd = 2;  // Start and end
    for (int i = 0; i < nk-1; i++) 
        tpd += (int) (pd[i] + gd);
    tpd += pd[nk-1];
    // -----------------------------------
    
    // Outgoing brepository ---------------
    
    rf   = Matrix<cxfl>  (tpd, nc);
    grad = Matrix<float> (tpd, 3);
    // -----------------------------------
    
    // RF Timing -------------------------

    for (int rc = 0; rc < nc; rc++) {
        
        int i = 1;
        
        for (int k = 0; k < nk; k++) {
            
            // RF action
            for (int p = 0; p < pd[k]; p++, i++) 
                rf (i,rc) = solution(k + rc*nk) / (float) pd[k] * cxfl(1000.0,0.0);

            // Gradient action, no RF
            if (k < nk-1)
                for (int g = 0; g <    gd; g++, i++)
                    grad (i,0) = 0.0;

        }
        
    }
    // -----------------------------------
    
    // Gradient and slew -----------------
    
    float gr = 0.0;
    float sr = 0.0;

    for (int gc = 0; gc < 3; gc++) {
        
        int i = 1;
        
        for (int k = 0; k < nk-1; k++) {
            
            // RF action, no gradients
            for (int p = 0; p < pd[k]; p++, i++) 
                grad (i,gc) = 0.0; 
            
            if (k < nk-1) 
                
                for (int g = 0; g <    gd; g++, i++) {
                    
                    sr = (k+1 < nk) ? ks(gc,k+1) - ks(gc,k) : - ks(gc,i);
                    sr = 4.0 * sr / (2 * PI * 4.2576e7 * gd * gd * 1.0e-5 * 1.0e-5);
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
    
    Matrix<cxfl> treg = lambda * eye<cxfl>(size(m,1));
    Matrix<cxfl> minv;
    
    minv  = m.prodt (m); 
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

        if ((gc && j && breakearly && res[gc] > res[gc-1]) || res[gc] < conv) {
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
        
        for (int j = 0; j < nc; j++)
            if (max_rf[i] < abs (solution[i+nk*j]) / (float)(10.0*pd[i])) 
                max_rf[i] = abs (solution[i+nk*j]) / (float)(10.0*pd[i]);

        max_rf[i] *= 100.0;

    }
    
	for (int i = 0; i < nk; i++)
		if (max_rf[i] > rflim) amps_ok = false;
	
	// Update Pulse durations if necessary -------
    
	if (!amps_ok) {
		
		printf ("Pulse amplitudes to high!\n  Updating pulse durations ... to "); fflush(stdout);
        
		for (int i = 0; i < nk; i++) {
			pd[i] = 1 + (int) (max_rf[i] * pd[i] / rflim); 
			printf ("%i ", 10*pd[i]); fflush(stdout);
		}
        
		printf ("[us] ... done.\n\n");
		
	} else 

		printf ("OK\n\n");

	return amps_ok;
	
}
