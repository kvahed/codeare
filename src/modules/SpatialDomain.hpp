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

#ifndef __SPATIAL_DOMAIN_HPP__
#define __SPATIAL_DOMAIN_HPP__

const int    GRAD_RASTER = 10;
const double GAMMA       = 42576000;

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

        int    m_nc;     /**<Transmit channels    */
        int*   m_pd;     /**< Pulse durations     */
        int    m_gd;     /**< Gradient durations  */
        int    m_ns;     /**< # Spatial positions */
        int    m_nk;     /**< # kt-points         */
        int    m_maxiter; /**< # Variable exchange method iterations */

        double m_lambda; /**< Tikhonov parameter  */
        double m_rflim;
        double m_conv;
        
        float* m_max_rf;

        std::string m_orient; /**< Orientation*/ 

    };

}
#endif /* __DUMMY_RECON_H__ */


/**
 * @brief           Normalised root-means-squared error
 *
 * @param  target   Target magnetisation
 * @param  result   Achieved result
 * @param           NRMSE
 */
void
NRMSE                         (const Matrix<raw>* target, const Matrix<raw>* result, float* nrmse) {

    float q = 0.0, n = 0.0;
    
    for (int i=0; i < target->Size(); i++)
        q += pow(abs(target->at(i)) - abs(result->at(i)), 2.0);
    
    q = sqrt(q)/target->norm().real();
    
    printf (" %.3f\n", q);
    
    nrmse[0] = 100.0 * q;

}


/**
 * @brief           Phase correction from off-resonance
 *
 * @param  target   Target magnetisation
 * @param  result   Achieved result
 * @return          Phase corrected 
 */
void
PhaseCorrection (Matrix<raw>* target, const Matrix<raw>* result) {
    
#pragma omp parallel default (shared) 
    {
        
        int tid      = omp_get_thread_num();
        int chunk    = target->Size() / omp_get_num_threads();
        
#pragma omp for schedule (dynamic, chunk)

        for (int i=0; i < target->Size(); i++) 
            if (abs(target->at(i)) > 0)
                target->at(i) = abs(target->at(i)) * result->at(i) / abs(result->at(i));
            else                
                target->at(i) = raw (0,0);    
        
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
void 
RFLimits            (const Matrix<raw>* solution, const int* pd, const int nk, const int nc, float* limits) {
    
    for (int i = 0; i < nk; i++) {

        limits[i] = 0.0;
        
        for (int j = 0; j < nc; j++)
            if (limits[i] < abs (solution->at(i+nk*j)) / pd[i]) 
                limits[i] = abs (solution->at(i+nk*j)) / pd[i];

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
 * @param  m        Out: m_xy
 */
void
STA (const Matrix<double>* ks, const Matrix<double>* r, const Matrix<raw>* b1, const Matrix<short>* b0, 
     const int             nc, const int            nk, const int          ns, const int            gd, Matrix<raw>* m) {
    
    raw         pgd = raw (0, 2.0 * PI * GAMMA * GRAD_RASTER); // 2* i * \pi * \gamma * \delta t

#pragma omp parallel default (shared) 
    {
        
        int tid      = omp_get_thread_num();
        int chunk    = nc / omp_get_num_threads();
        
#pragma omp for schedule (dynamic, chunk)
        
        for (int c = 0; c < nc; c++) 
            for (int k = 0; k < nk; k++) 
                for (int s = 0; s < ns; s++) 
                    m->at (c*nk*ns + k*ns + s) = 
                        pgd * b1->at(c*ns + s) *                           // b1 (s,c)
                        exp (raw(0, 2.0 * PI * gd * (float) b0->at(s))) *  // off resonance: exp (2i\pidb0dt)  
                        exp (raw(0,(ks->at(k)*r->at(s) + ks->at(k+nk)*r->at(s+ns) + ks->at(k+2*nk)*r->at(s+2*ns)))); // encoding: exp (i k(t) r)
        
    }
    
}


