#include "CPUSimulator.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "linalg/Lapack.hpp"
#include "mri/MRI.hpp"

using namespace RRStrategy;


/*
#include "adolc.h"

inline static void 
DiffRot (const Matrix<float> n, Matrix<float>& lm) {

	adouble n1, n2, n3;
	adouble in1, in2, in3;
	adouble out1, out2, out3;

}
*/


/**
 * @brief       Rotate magnetisation around rotation axis
 *
 * @param  n    Rotation axis
 * @param  lm   On entry magnetisation vector to be rotated.<br/>
 *              On exit rotated magnetisation vector.
 */
inline static void 
Rotate (const Matrix<float>& n, Matrix<float>& lm) {
    
	float r  [9];
	float nm [3];

    float phi = sqrt(n[X]*n[X] + n[Y]*n[Y] + n[Z]*n[Z]);

    // Identity
    if (!phi) {

		r[0] = 1.0; r[3] = 0.0; r[6] = 0.0;
		r[1] = 0.0; r[4] = 1.0; r[7] = 0.0;
        r[2] = 0.0; r[5] = 0.0; r[8] = 1.0;
    
	} else {

        // Cayley-Klein parameters
        float hp    =  0.5    *phi;        
        float sp    =  sin(hp)/phi; /* /phi because n is unit length in defs. */
        float ar    =  cos(hp);
        float ai    = -n[Z]*sp;
        float br    =  n[Y]*sp;
        float bi    = -n[X]*sp;
    
        float arar  = ar*ar;
        float aiai  = ai*ai;
        float arai2 = ar*ai*2.0;
        float brbr  = br*br;
        float bibi  = bi*bi;
        float brbi2 = br*bi*2.0;
        float arbi2 = ar*bi*2.0;
        float aibr2 = ai*br*2.0;
        float arbr2 = ar*br*2.0;
        float aibi2 = ai*bi*2.0;
        
		float h1    = arar - aiai;
		float h2    = bibi - brbr;

        r[0] =  h1    + h2;
        r[1] = -arai2 - brbi2;
        r[2] = -arbr2 + aibi2;
        r[3] =  arai2 - brbi2; 
        r[4] =  h1    - h2;
        r[5] = -aibr2 - arbi2;
        r[6] =  arbr2 + aibi2;
        r[7] =  arbi2 - aibr2;
        r[8] =  arar  + aiai - brbr - bibi;
        
    }

	nm[X] = r[0]*lm[X] + r[3]*lm[Y] + r[6]*lm[Z];
	nm[Y] = r[1]*lm[X] + r[4]*lm[Y] + r[7]*lm[Z];
	nm[Z] = r[2]*lm[X] + r[5]*lm[Y] + r[8]*lm[Z];

	lm[0] = nm[0];
	lm[1] = nm[1];
	lm[2] = nm[2];

}



void 
SimulateAcq (const Matrix<cxfl>&  b1, const Matrix<float>&  gr, const Matrix<float>&   r, 
             const Matrix<float>& b0, const Matrix<cxfl>&  mt0, const Matrix<float>& ml0,
			 const Matrix<float>& ic,
			 const int&           np, const float&          dt, const bool&            v, 
			 const size_t&        nc, const size_t&         nt, const float&         gdt,        
			       Matrix<cxfl>&  rf) {


	size_t            nr   = r.Dim(1);
	float             nrs  = (float) nr;

    ticks             tic  = getticks();

	Matrix<cxfl>      sig (nt,nc,np); /*<! Signal repository  */
	
#pragma omp parallel
    {
		
		Matrix<float> n   ( 3,1);  // Rotation axis
		Matrix<float> lm  ( 3,1);  // Magnetisation
		Matrix<float> tmp ( 3,1);  // Temporary magnetisation
		
		Matrix<float> lr  ( 3,1);  // Local spatial vector
		Matrix<cxfl>  ls  (nc,1);  // Local sensitivity
		float         lb0;
		
        omp_set_num_threads(np);
		
#pragma omp for schedule (guided) 
		
		for (size_t pos = 0; pos < nr; pos++) {

			lm[X] = mt0[pos].real()*ic[pos]; 
			lm[Y] = mt0[pos].imag()*ic[pos]; 
			lm[Z] = ml0[pos]       *ic[pos];

			if ((lm[X]+lm[Y]+lm[Z]) > 0.0) {

				lb0   = b0[pos]*TWOPI;

				lr[0] = r[pos*3];
				lr[1] = r[pos*3+1];
				lr[2] = r[pos*3+2];

				for (size_t c = 0; c < nc; c++) 
					ls[c] = conj(b1(pos,c));
				
				// Run over time points
				size_t t = nt;
				while (t--) {
					
					tmp[0] = lm[0];
					tmp[1] = lm[1];
					tmp[2] = lm[2];
					
					// Rotate axis (only gradients)
					n[2] = - gdt * (gr(X,t)*lr[X] + gr(Y,t)*lr[Y] + gr(Z,t)*lr[Z] - lb0);
					
					Rotate (n, lm);
					
					// Weighted contribution to all coils
					cxfl mxy ((lm[X]+tmp[X])/2,(lm[Y]+tmp[Y])/2); 
					for (size_t c = 0; c < nc; c++)
						sig.At(t,c,omp_get_thread_num()) += ls[c]*mxy;
					
				}
				
			}
			
		}
		
#pragma omp for  schedule (guided) 
		for (size_t i = 0; i < nt*nc; i++) {
			rf[i] = cxfl(0.0,0.0);
			for (int p = 0; p < np; p++)
				rf[i] += sig[p*nt*nc+i];
			rf[i] /= nrs;
		}
		
	}
	
	if (v) printf ("(a: %.4fs)", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate()); fflush(stdout);
	
}



void 
SimulateExc  (const Matrix<cxfl>&   b1, const Matrix<float>&  gr, const Matrix< cxfl>& rf, 
              const Matrix<float>&   r, const Matrix<float>&  b0, const Matrix<cxfl>& mt0, 
			  const Matrix<float>& ml0, const Matrix<float>& jac, 
			  const size_t&         np, const float&          dt, const bool&           v, 
			  const size_t&         nc, const size_t&         nt, const float&        gdt, 
			        Matrix<cxfl>&  mxy,       Matrix<float>&  mz) {
    
	size_t            nr   = r.Dim(1);
	
    ticks             tic  = getticks();

#pragma omp parallel default(shared)
	{
		
		Matrix<float> n   ( 3,1);  // Rotation axis
		Matrix<float> lm  ( 3,1);  // Magnetisation
		Matrix<float> lr  ( 3,1);  // Local spatial vector
		Matrix<cxfl>  ls  (nc,1);  // Local sensitivity
		float         lb0;
		
		omp_set_num_threads(np);
	
#pragma omp for schedule (guided) 
	
		for (size_t pos = 0; pos < nr; pos++) {
			
			// Start with equilibrium
			lm[X] = mt0[pos].real();
			lm[Y] = mt0[pos].imag();
			lm[Z] = ml0[pos];
			
			if ((lm[X]+lm[Y]+lm[Z]) > 0.0) {

				lb0 = b0[pos]*TWOPI;

				lr[0] = r[pos*3  ];
				lr[1] = r[pos*3+1];
				lr[2] = r[pos*3+2];

				for (size_t i = 0; i < nc; i++) ls[i] = b1 (pos,i);
				
				// Time points
				for (size_t t = 0; t < nt; t++) {
					
					cxfl rfs = cxfl (0.0,0.0);
					for (size_t i = 0; i < nc; i++) 
						rfs += rf(t,i)*ls[i];
					rfs *= jac[t];
					
					n[0]  = gdt * -rfs.imag() * 1.0e-3;
					n[1]  = gdt *  rfs.real() * 1.0e-3;
					n[2]  = gdt * (gr(X,t) * lr[X] + gr(Y,t) * lr[Y] + gr(Z,t) * lr[Z] - lb0) ;
					
					Rotate (n, lm);
					
				}
				
				mxy [pos] = cxfl (lm[X], lm[Y]);
				mz  [pos] = lm[Z]; 
				
			}
			
		}
		
	}
	
	if (v) printf ("(e: %.4fs)", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate()); fflush(stdout);
	
}



CPUSimulator::CPUSimulator (SimulationBundle* sb) {
    
    m_sb  = sb;

    m_gdt  = GAMMA * TWOPI * m_sb->dt;
    m_nt   = m_sb->g->Dim(1);    // Time points
    m_nc   = m_sb->b1->Dim(1);     // # channels
	m_nr   = m_sb->r->Dim(1);      // # spatial positions
	m_rfsc = 1.0;
	
	// Intensity map for correction
	m_ic = IntensityMap (*(m_sb->b1), false);

	if (m_sb->roi->Size() == 1)
		m_sb->roi = m_sb->smz;

}



CPUSimulator::~CPUSimulator () {}


void 
CPUSimulator::Simulate () {

	Matrix<cxfl>&   b1 = *(m_sb->b1);
	Matrix<float>&   g = *(m_sb->g);
	Matrix<float>&  rv = *(m_sb->r);
	Matrix<float>&  b0 = *(m_sb->b0);
	Matrix<cxfl>& tmxy = *(m_sb->tmxy);
	Matrix<float>& tmz = *(m_sb->tmz);
	Matrix<cxfl>&   rf = *(m_sb->rf);
	Matrix<cxfl>& smxy = *(m_sb->smxy);
	Matrix<float>& smz = *(m_sb->smz);
	Matrix<float>& roi = *(m_sb->roi);
	Matrix<float>& jac = *(m_sb->jac);
	Matrix<cxfl>&  mxy = *(m_sb->mxy);
	Matrix<float>&  mz = *(m_sb->mz);

	int             np = m_sb->np;
	float           dt = m_sb->dt;
	float           tl = m_sb->lambda;
	vector<float>  res;	

	bool           cb0 = m_sb->cb0;

	SimulateAcq (b1, g, rv, (cb0) ? b0 : zeros<float> (m_nr,1), tmxy, tmz, m_ic, np, dt, false, m_nc, m_nt, m_gdt, rf);

	if (m_sb->mode) {                                   // Optimise

		int   iters = 0;
		float rn = 0.0, rno = 0.0, an = 0.0;
		cxfl  rtmp = cxfl(0.0,0.0);
		Matrix<cxfl> p, r, q, a;

		p = rf; q = p; r = p;
		a = zeros<cxfl> (size(rf));

		an = pow(norm(p), 2.0);
		rn = an;
		
		for (iters = 0; iters < m_sb->cgit; iters++) { 	// CG loop
			
			res.push_back(rn/an);

			if (iters == 0) printf ("                        "); 
			printf ("  %03i: CG residuum: %.9f\n", iters, res.at(iters));
			if (res.at(iters) <= m_sb->cgeps) break;    // Convergence?

			SimulateExc (b1, g, p, rv, b0, smxy, roi,  jac, np, dt, true, m_nc, m_nt, m_gdt, mxy, mz); // E^H
			SimulateAcq (b1, g,    rv, b0,  mxy,  mz, m_ic, np, dt, true, m_nc, m_nt, m_gdt,       q); // E
			q += tl * p; // L_2 penalty on solution

			rtmp = (rn / p.dotc(q));
			
			a   += rtmp * p;
			r   -= rtmp * q;
			rno  = rn;
			rn   = pow(norm(r),2.0);
			p   *= rn/rno;
			p   += r;
			
		}

		rf = a; // Final pulses
		SimulateExc (b1, g, rf, rv, b0, smxy, smz, jac, np, dt, false, m_nc, m_nt, m_gdt, mxy, mz); // Excitation 
		
	}

}


