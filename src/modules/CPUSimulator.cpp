#include "CPUSimulator.hpp"

using namespace RRStrategy;

/**
 * @brief       Rotate magnetisation around rotation axis
 *
 * @param  n    Rotation axis
 * @param  r    Magnetisation vector to be rotated
 */
inline void 
Rotate (const Matrix<float>& n, Matrix<float>& lm) {
    
	float r   [9];
	float tmp [3];

    float phi = sqrt(n[X]*n[X] + n[Y]*n[Y] + n[Z]*n[Z]);

    // Identity
    if (!phi) {

		r[0] = 1.0;
        r[1] = 0.0;
        r[2] = 0.0;
        r[3] = 0.0;
        r[4] = 1.0;
        r[5] = 0.0;
        r[6] = 0.0;
        r[7] = 0.0;
        r[8] = 1.0;
    
	} else {

        // Cayley-Klein parameters
        float hp    =  phi/2;        
        float cp    =  cos(hp);
        float sp    =  sin(hp)/phi; /* /phi because n is unit length in defs. */
        float ar    =  cp;
        float ai    = -n[Z]*sp;
        float br    =  n[Y]*sp;
        float bi    = -n[X]*sp;
    
        float arar  =   ar*ar;
        float aiai  =   ai*ai;
        float arai2 = 2*ar*ai;
        float brbr  =   br*br;
        float bibi  =   bi*bi;
        float brbi2 = 2*br*bi;
        float arbi2 = 2*ar*bi;
        float aibr2 = 2*ai*br;
        float arbr2 = 2*ar*br;
        float aibi2 = 2*ai*bi;
        
        r[0] =  arar  - aiai - brbr + bibi;
        r[1] = -arai2 - brbi2;
        r[2] = -arbr2 + aibi2;
        r[3] =  arai2 - brbi2; 
        r[4] =  arar  - aiai + brbr - bibi;
        r[5] = -aibr2 - arbi2;
        r[6] =  arbr2 + aibi2;
        r[7] =  arbi2 - aibr2;
        r[8] =  arar  + aiai - brbr - bibi;
        
    }

	tmp[X] = r[0]*lm[X] + r[3]*lm[Y] + r[6]*lm[Z];
	tmp[Y] = r[1]*lm[X] + r[4]*lm[Y] + r[7]*lm[Z];
	tmp[Z] = r[2]*lm[X] + r[5]*lm[Y] + r[8]*lm[Z];

	lm[0] = tmp[0];
	lm[1] = tmp[1];
	lm[2] = tmp[2];

}



void 
SimulateAcq (const Matrix<cplx>&  b1, const Matrix<float>&  gr, const Matrix<float>&   r, 
             const Matrix<float>& b0, const Matrix<cplx>&  mt0, const Matrix<float>& ml0,
			 const int&           np, const float&          dt, const bool&            v, 
			 const size_t&        nc, const size_t&         nt, const float&         gdt,        
			       Matrix<cplx>&  rf) {


	size_t            nr   = r.Dim(1);

    ticks             tic  = getticks();
	printf ("    %02i-ch acquisition of %i isochromats ... ", nc, nr); fflush(stdout);

	Matrix<cplx>   sig (nt,nc,np); /*<! Signal repository  */
	
#pragma omp parallel
    {
		
		Matrix<float> n   ( 3,1);  // Rotation axis
		Matrix<float> lm  ( 3,1);  // Magnetisation
		Matrix<float> tmp ( 3,1);  // Temporary magnetisation
		
		Matrix<float> lr  ( 3,1);  // Local spatial vector
		Matrix<cplx>  ls  (nc,1);  // Local sensitivity
		
        omp_set_num_threads(np);
		
#pragma omp for schedule (guided) 
		
		for (size_t pos = 0; pos < nr; pos++) {

			lm[X] = mt0[pos].real(); 
			lm[Y] = mt0[pos].imag(); 
			lm[Z] = ml0[pos];
			
			if ((lm[X]+lm[Y]+lm[Z]) > 0.0) {

				for (size_t c = 0; c < nc; c++) ls[c] = conj(b1(pos,c));
				for (size_t i = 0; i <  3; i++) lr[i] =       r(i,pos);
				
				// Run over time points
				for (size_t t = 0; t < nt; t++) {
					
					tmp[0] = lm[0];
					tmp[1] = lm[1];
					tmp[2] = lm[2];
					
					// Rotate axis (only gradients)
					n[2] = - gdt * (gr(X,t)*lr[X] + gr(Y,t)*lr[Y] + gr(Z,t)*lr[Z] + b0[pos]*TWOPI);
					
					Rotate (n, lm);
					
					// Weighted contribution to all coils
					cplx mxy  = cplx((lm[X]+tmp[X])/2,(lm[Y]+tmp[Y])/2); 
					cplx scon = 0.0;
					
					for (size_t c = 0; c < nc; c++)
						sig.At(t,c,omp_get_thread_num()) += ls[c]*mxy;
					
				}
				
			}
			
		}

#pragma omp for  schedule (guided) 
		for (size_t i = 0; i < nt*nc; i++) {
			rf[i] = cplx(0.0,0.0);
			for (int p = 0; p < np; p++)
				rf[i] += sig[p*nt*nc+i];
			rf[i] /= nr;
			//rf[i] *= 37.3814;
		}
		
	}

	printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

}



void 
SimulateExc  (const Matrix<cplx>&   b1, const Matrix<float>&  gr, const Matrix< cplx>& rf, 
              const Matrix<float>&   r, const Matrix<float>&  b0, const Matrix<cplx>& mt0, 
			  const Matrix<float>& ml0, const Matrix<float>& jac, const size_t&        np, 
			  const float&          dt, const bool&            v, const size_t&        nc, 
			  const size_t&         nt, const float&         gdt,       
			        Matrix<cplx>&  mxy,       Matrix<float>&  mz) {
    
	size_t            nr   = r.Dim(1);
	
    ticks             tic  = getticks();

	printf ("    %02i-ch excitation of %i isochromats ... ", nc, nr); fflush(stdout);

#pragma omp parallel default(shared)
	{
		
		Matrix<float> n   ( 3,1);  // Rotation axis
		Matrix<float> lm  ( 3,1);  // Magnetisation
		Matrix<float> lr  ( 3,1);  // Local spatial vector
		Matrix<cplx>  ls  (nc,1);  // Local sensitivity
		
		omp_set_num_threads(np);
	
#pragma omp for schedule (guided) 
	
		for (size_t pos = 0; pos < nr; pos++) {
			
			// Start with equilibrium
			lm[X] = mt0[pos].real();
			lm[Y] = mt0[pos].imag();
			lm[Z] = ml0[pos];
			
			if ((lm[X]+lm[Y]+lm[Z]) > 0.0) {
				
				for (size_t i = 0; i <  3; i++) lr[i] = r  (i,pos);
				for (size_t i = 0; i < nc; i++) ls[i] = b1 (pos,i);
				
				// Time points
				for (size_t t = 0; t < nt; t++) {
					
					size_t rt = nt-1-t;
					
					cplx rfs = cplx (0.0,0.0);
					for (size_t i = 0; i < nc; i++) 
						rfs += rf(rt,i)*ls[i];
					rfs *= jac[rt];
					
					n[0] = gdt * -rfs.imag();
					n[1] = gdt *  rfs.real();
					n[2] = gdt * (gr(X,rt) * lr[X] + gr(Y,rt) * lr[Y] + gr(Z,rt) * lr[Z] - b0[pos]*TWOPI) ;
					
					Rotate (n, lm);
					
				}
				
				mxy [pos] = cplx (lm[X], lm[Y]);
				mz  [pos] = lm[Z]; 
				
			}
			
		}
		
	}
	
	printf ("done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

}



CPUSimulator::CPUSimulator (SimulationBundle* sb) {
    
    m_sb  = sb;

    m_gdt = GAMMARAD * m_sb->dt;
    m_nt  = m_sb->agr->Dim(1);      // Time points
    m_nc  = m_sb->b1->Dim(1);     // # channels
	
	m_sb->Dump("sb.mat");

}



CPUSimulator::~CPUSimulator () {}


void 
CPUSimulator::Simulate () {

    ticks            tic  = getticks();                                 // Start timing
	
	SimulateAcq (*(m_sb->b1), *(m_sb->agr), *(m_sb->r), *(m_sb->b0), 
				 *(m_sb->tmxy), *(m_sb->tmz), m_sb->np, m_sb->dt, m_sb->v, 
				 m_nc, m_nt, m_gdt, *(m_sb->rf));
	if (!m_sb->exc) return;

	SimulateExc (*(m_sb->b1), *(m_sb->agr), *(m_sb->rf), *(m_sb->r), 
				 *(m_sb->b0), *(m_sb->smxy), *(m_sb->smz), *(m_sb->jac), 
				 m_sb->np, m_sb->dt, m_sb->v, m_nc, m_nt, m_gdt, 
				 *(m_sb->mxy), *(m_sb->mz));

	if (m_sb->mode) {

		int   iters = 0;
		
		float rn    = 0.0;
		float an    = 0.0;
		cplx  rtmp  = cplx(0.0,0.0);
		
		Matrix<cplx> p, r, q, a;

		p = *(m_sb->mxy);
		q = p;
		r = p;
		
		vector<float> res;
		
		an = pow(p.Norm().real(), 2);
		
		char buffer [7];
		std::string s;
		// CGNR loop
		for (iters = 0; iters < m_sb->cgit; iters++) {
			
			rn = pow(r.Norm().real(), 2);
			res.push_back(rn/an);
			printf ("  %03i: CG residuum: %.9f\n", iters, res.at(iters));
			
			// Convergence ? ----------------------------------------------
			if (res.at(iters) <= m_sb->cgeps) break;

			sprintf(buffer, "q%02i.mat", iters);
			s = std::string (buffer);
			q.MXDump(s.c_str(), "q");

			SimulateAcq (*(m_sb->b1), *(m_sb->agr), *(m_sb->r), *(m_sb->b0), 
						 q, *(m_sb->tmz), m_sb->np, m_sb->dt, m_sb->v, 
						 m_nc, m_nt, m_gdt, *(m_sb->rf));
			
			SimulateExc (*(m_sb->b1), *(m_sb->agr), *(m_sb->rf), *(m_sb->r), 
						 *(m_sb->b0), *(m_sb->smxy), *(m_sb->smz), *(m_sb->jac), 
						 m_sb->np, m_sb->dt, m_sb->v, m_nc, m_nt, m_gdt, 
						 p, *(m_sb->mz));

			rtmp  = (rn / (p.dotc(q)));
			if (iters == 0)
				a     = (p    * rtmp);
			else
				a    += (p    * rtmp);
			r    -= (q    * rtmp);
			p    *= cplx(pow(r.Norm().real(), 2)/rn);
			p    += r;
			
			sprintf(buffer, "a%02i.mat", iters);
			s  = std::string(buffer);
			a.MXDump(s.c_str(), "a");

		}

	}

}


