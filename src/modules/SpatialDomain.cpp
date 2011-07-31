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

const int    GRAD_RASTER = 10;
const double GAMMA       = 42576000;


using namespace RRStrategy;


/**
 * @brief           Phase correction from off-resonance
 *
 * @param  target   Target magnetisation
 * @param  result   Achieved result
 * @return          Phase corrected 
 */
void
PhaseCorrection     (Matrix<raw>* target, const Matrix<raw>* result);

/**
 * @brief           Normalised root-means-squared error
 *
 * @param  target   Target magnetisation
 * @param  result   Achieved result
 * @return          NRMSE
 */
float       
NRMSE               (const Matrix<raw>* target, const Matrix<raw>* result);

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
RFLimits            (const Matrix<raw>* solution, const int* pd, const int nk, const int nc, float* limits);

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
 * @return          magnetisation
 */
Matrix<raw>
STA                 (const Matrix<double>* ks, const Matrix<double>* r, const Matrix<raw>* b1, const Matrix<short>* b0, 
					 const int             nc, const int            nk, const int          ns, const int            gd);




SpatialDomain::SpatialDomain  () {}


SpatialDomain::~SpatialDomain () {}


RRSModule::error_code
SpatialDomain::Init           ()  {

	printf ("Intialising SpatialDomain ...\n");

    RRSModule::error_code e = OK;

	// # kt points ---------------------------
    Attribute ("Nk",      &m_nk);
	printf ("  # kt-points: %i \n", m_nk);

    m_max_rf = (float*) malloc (m_nk * sizeof(float));
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
    // m_raw:    Pulses
    // m_helper: Pulse durations
    // ----------------------------
    Matrix<raw> treg         = Matrix<raw>::id(m_nc * m_nk) * raw (m_lambda, 0);

    bool        pulse_amp_ok = false;

    Matrix<raw> solution;
    Matrix<raw> tmp;
    Matrix<raw> final;    

    // Start clock ------------------------
    ticks vestart = getticks();
	
	while (!pulse_amp_ok) {
		
        Matrix<raw> m    = STA (&m_kspace, &m_helper, &m_rhelper, &m_pixel, m_nc, m_nk, m_ns, m_gd);
		Matrix<raw> minv = m.tr();

		minv  = minv.prod (m);
        minv  = minv + treg;
		minv  = minv.Pinv();
        minv  = minv.prodt (m);
        
        // Valriable exchange method --------------

        std::vector<double> res;
        
        printf ("  Starting variable exchange method ...\n");

        for (int j = 0; j < m_maxiter; j++) {

            solution = minv.prod(m_raw);
            tmp      = m.prod(solution);

			printf ("  %03i: ", j);
            res.push_back (NRMSE (&m_raw, &tmp));
			
			if (res.at(j) < m_conv) break;
            
            final    = solution;
			PhaseCorrection (&m_raw, &tmp);
			
		}
		
		printf ("... done. Checking pulse amplitudes ... \n");
		
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
		
	}

	printf ("... done. WTime: %.4f seconds.\n", elapsed(getticks(), vestart) / ClockRate());

	m_raw = final;

    return RRSModule::OK;

}



float      
NRMSE                         (const Matrix<raw>* target, const Matrix<raw>* result) {

	float q = 0.0, n = 0.0;
	
	for (int i=0; i < target->Size(); i++)
		q += pow(abs(target->at(i)) - abs(result->at(i)), 2.0);
	
	q = sqrt(q)/target->norm().real();
	
    printf (" %.3f\n", q);
	
	return 100.0 * q;

}



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


void 
RFLimits            (const Matrix<raw>* solution, const int* pd, const int nk, const int nc, float* limits) {
    
    for (int i = 0; i < nk; i++) {

        limits[i] = 0.0;
        
        for (int j = 0; j < nc; j++)
            if (limits[i] < abs (solution->at(i+nk*j)) / pd[i]) 
                limits[i] = abs (solution->at(i+nk*j)) / pd[i];

    }
        
}


Matrix<raw>
STA (const Matrix<double>* ks, const Matrix<double>* r, const Matrix<raw>* b1, const Matrix<short>* b0, 
	 const int             nc, const int            nk, const int          ns, const int            gd) {
    
    Matrix<raw> mxy (ns, nk*nc);
    raw         pgd = raw (0, 2.0 * PI * GAMMA * GRAD_RASTER);            // 2* i * \pi * \gamma * \delta t

#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		int chunk    = nc / omp_get_num_threads();
		
#pragma omp for schedule (dynamic, chunk)
		
		for (int c = 0; c < nc; c++) 
			for (int k = 0; k < nk; k++) 
				for (int s = 0; s < ns; s++) 
					mxy.at (c*nk*ns + k*ns + s) = 
						pgd * b1->at(c*ns + s) *                           // b1 (s,c)
						exp (raw(0, 2.0 * PI * gd * (float) b0->at(s))) *  // off resonance: exp (2i\pidb0dt)  
						exp (raw(0,(ks->at(k)*r->at(s) + ks->at(k+nk)*r->at(s+ns) + ks->at(k+2*nk)*r->at(s+2*ns)))); // encoding: exp (i k(t) r)

	}

	return mxy;
	
}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {

    return new SpatialDomain;

}



extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {

    delete p;

}
