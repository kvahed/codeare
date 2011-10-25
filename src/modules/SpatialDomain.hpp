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

#ifndef __SPATIAL_DOMAIN_HPP__
#define __SPATIAL_DOMAIN_HPP__

#include "ReconStrategy.hpp"

using namespace RRServer;


/**
 * @brief Reconstruction startegies
 */
namespace RRStrategy {

    /**
     * @brief Spatial domain method for RF pulse generation
     */
    class SpatialDomain : public ReconStrategy {
        
        
    public:
        
        /**
         * @brief Default constructor
         */
        SpatialDomain  ();
        
        /**
         * @brief Default destructor
         */
        virtual 
        ~SpatialDomain ();
        
        
        /**
         * @brief Process
         */
        virtual RRSModule::error_code
        Process ();

        
        /**
         * @brief Initialise
         */
        virtual RRSModule::error_code
        Init ();
        

        /**
         * @brief Finalise
         */
        virtual RRSModule::error_code
        Finalise ();



    private:
 
        int         m_nc;       /**< Transmit channels       */
        int*        m_pd;       /**< Pulse durations         */
        int         m_gd;       /**< Gradient durations      */
        int         m_ns;       /**< # Spatial positions     */
        int         m_nk;       /**< # kt-points             */
        int         m_maxiter;  /**< # Variable exchange method iterations */
		int         m_verbose;  /**< Verbose output. All intermediate results. */
		int         m_breakearly;  /**< Verbose output. All intermediate results. */

        double      m_lambda;   /**< Tikhonov parameter      */
        double      m_rflim;    /**< Maximum rf amplitude    */
        double      m_conv;     /**< Convergence criterium   */
        
        float*      m_max_rf;   /**< Maximum reached RF amps */

        std::string m_orient;   /**< Orientation             */ 
		std::string m_ptxfname; /**< PTX file name           */

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
inline void
NRMSE                         (const Matrix<cplx>& target, const Matrix<cplx>& result, const int& iter, float& nrmse) {

    float q = 0.0;
    
    for (int i=0; i < target.Size(); i++)
        q += pow(abs(target[i]) - abs(result[i]), 2);
    
    q = sqrt(q)/target.Norm().real();
    
	if (iter % 5 == 0 && iter > 0)
		printf ("\n");
    printf ("    %03i %.6f", iter, q);

    nrmse = 100.0 * q;

}


/**
 * @brief           Phase correction from off-resonance
 *
 * @param  target   Target magnetisation
 * @param  result   Achieved result
 * @return          Phase corrected 
 */
inline void
PhaseCorrection (Matrix<cplx>& target, const Matrix<cplx>& result) {
    
	size_t n = target.Size();

#pragma omp parallel default (shared) 
    {
        
        int tid      = omp_get_thread_num();
        int chunk    = n / omp_get_num_threads();
        
#pragma omp for schedule (dynamic, chunk)

        for (size_t i = 0; i < n; i++) 
            target[i] = (abs(target[i]) > 0) ? abs(target[i]) * result[i] / abs(result[i]) :  cplx(0,0);    
        
    }
    
}


/**
 * @brief           RF limts
 *
 * @param  solution In:  Calculated solution
 * @param  pd       In:  Pulse durations
 * @param  nk       In:  # Pulses
 * @param  nc       In:  # Coils
 * @param  limits   Out: limits
 */
inline void 
RFLimits            (const Matrix<cplx>& solution, const int* pd, const int& nk, const int& nc, float* limits) {
    
    for (int i = 0; i < nk; i++) {

        limits[i] = 0.0;
        
        for (int j = 0; j < nc; j++)
            if (limits[i] < abs (solution[i+nk*j]) / (float)(10.0*pd[i])) 
                limits[i] = abs (solution[i+nk*j]) / (float)(10.0*pd[i]);

		limits[i] *= 100.0;

    }

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
void
STA (const Matrix<double>& ks, const Matrix<double>& r, const Matrix<cplx>& b1, const Matrix<short>& b0, 
	 const int& nc, const int& nk, const int& ns, const int gd, const int* pd, Matrix<cplx>& m) {
    
	float* d = (float*) malloc (nk * sizeof (float));
	float* t = (float*) malloc (nk * sizeof (float));

	for (int i = 0; i< nk; i++)
		d[i] = (i==0) ? pd[i] + gd : d[i-1] + pd[i] + gd;

	for (int i = 0; i< nk; i++)		
		t[i] = d [nk-i-1];
	
	for (int i = 0; i < nk-1; i++)
		d[i] = 1.0e-5 * t[i+1] + 1.0e-5 * pd[i]/2;

	d[nk-1] = 1.0e-5 * pd[nk-1] / 2;

    cplx pgd = cplx (0, 2.0 * PI * 4.2576e7 * 1.0e-5); 

#pragma omp parallel default (shared) 
    {
        
        int tid      = omp_get_thread_num();
        int chunk    = nc / omp_get_num_threads();
        
#pragma omp for// schedule (dynamic, chunk)
        
		// pTX STA 
        for (int c = 0; c < nc; c++) 
            for (int k = 0; k < nk; k++) 
                for (int s = 0; s < ns; s++) 
                    m(s, c*nk + k) = 
						// b1 (s,c)
                        pgd * b1(s,c) *
						// off resonance: exp (2i\pidb0dt)  
                        exp (cplx(0, 2.0 * PI * d[k] * (float) b0(s))) *
						// encoding: exp (i k(t) r)
                        exp (cplx(0,(ks(0,k)*r(0,s) + ks(1,k)*r(1,s) + ks(2,k)*r(2,s))));
        
    }
	
	free (t);
	free (d);
    
}


/**
 * @brief Construct actual pulses
 *
 *
 */
inline void
PTXTiming (const Matrix<cplx>& rf, const Matrix<double>& ks, const int* pd, const int& gd, const int& nk, const int& nc, Matrix<cplx>& timing) {

	// Total pulse duration --------------

	int tpd = 2;  // Start and end
	// Total excitation duration
	for (int i = 0; i < nk-1; i++) 
		tpd += (int) (pd[i] + gd);
	tpd += pd[nk-1];
	// -----------------------------------

	// Outgoing brepository ---------------

	// Time scale
	timing.Dim(COL) = tpd;    
	// Nc Channels + 3 gradients
	timing.Dim(LIN) = nc + 3; 
	timing.Reset();
	// -----------------------------------
	
	// RF Timing -------------------------

	for (int rc = 0; rc < nc; rc++) {
		
		int i = 1;
		
		for (int k = 0; k < nk; k++) {
			
			// RF action
			for (int p = 0; p < pd[k]; p++, i++) 
				//timing(i,rc) = conj(rf(k + rc*nk)) / (float)pd[k] * cplx(1000.0,0.0);
				timing(i,rc) = rf(k + rc*nk) / (float)pd[k] * cplx(1000.0,0.0);

			// Gradient action, no RF
			if (k < nk-1)
				for (int g = 0; g <    gd; g++, i++)
					timing(i,rc) = cplx (0.0, 0.0);

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
				timing(i,nc+gc) = 0.0; 

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
				
				timing(i,nc+gc) = gr; 
				
			}
			
		} 
		
	}
	// ----------------------------------

}
