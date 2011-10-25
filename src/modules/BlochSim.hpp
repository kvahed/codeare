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

#include "OMP.hpp"

#define GAMMARAD 2.6753e3
#define TWOPI	 6.283185

enum coords {
	X, Y, Z
};


/**
 * @brief       Rotation matrix around vector n
 *
 * @param  n    Rotation axis
 * @param  r    Rotation matrix (output)
 */
void RotationMatrix (const Matrix<double>& n, Matrix<double>& r) {
	
	double phi = sqrt(n[X]*n[X] + n[Y]*n[Y] + n[Z]*n[Z]);

	// Identity
	if (!phi)
		r = Matrix<double>::Id(3);
	
	else {

		// Cayley-Klein parameters
		double hp    =  phi/2;		
		double cp    =  cos(hp);
		double sp    =  sin(hp)/phi; /* /phi because n is unit length in defs. */
		double ar    =  cp;
		double ai    = -n[Z]*sp;
		double br    =  n[Y]*sp;
		double bi    = -n[X]*sp;
	
		double arar  =   ar*ar;
		double aiai  =   ai*ai;
		double arai2 = 2*ar*ai;
		double brbr  =   br*br;
		double bibi  =   bi*bi;
		double brbi2 = 2*br*bi;
		double arbi2 = 2*ar*bi;
		double aibr2 = 2*ai*br;
		double arbr2 = 2*ar*br;
		double aibi2 = 2*ai*bi;
		
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

}




/**
 * @brief       Simulate single shot excitation of a single isochromat of r matrix along gradient trajectory<br/>
 *              (i.e. inverse Fourier transform incl. effect of Receive and b0 maps)</br>
 *              Expects res to carry the correct size and dimensions
 *
 * @param  txm  Receive sensitivities
 * @param  rf   Complex RF field
 * @param  gr   Gradient trajectory
 * @param  r    Positions
 * @param  b0m  B0 map
 * @param  v    Verbose
 * @param  n    Position of isochromat in r
 * @param  m    Resulting magntisation (output)
 */
void SimulateExc  (const Matrix<cplx>&   txm, const Matrix<cplx>&  rf, const Matrix<double>& gr, const Matrix<double>& r, 
				   const Matrix<double>& b0m, const double&        dt, const bool&            v, const size_t&       pos, 
				   const int&            tid,       Matrix<double>& m) {

	
	assert (gr.Dim(1) == rf.Dim(0));
	
	size_t nt = gr.Dim(1);  // Time points
	size_t nc = txm.Dim(1); // # channels
	
	Matrix<double> n   ( 3,1);  // Rotation axis
	Matrix<double> rot ( 3,3);  // Rotation matrix
	Matrix<double> ml  ( 3,1);  // Local magnetisation
	Matrix<double> tmp ( 3,1);  // Local magnetisation

	Matrix<double> lr  (3,1);    // Local spatial vector
	Matrix<cplx>   ls  (nc,1);   // Local sensitivity

	for (size_t i = 0; i <  3; i++)
		lr[i] = r   (i,pos);
	for (size_t i = 0; i < nc; i++)
		ls[i] = txm(pos,i);

	ml[Z] = 1.0;

	double gdt = GAMMARAD * dt;
	
	// Run over time points
	for (size_t t = 0; t < nt; t++) {
		
		cplx rfs = cplx (0.0,0.0);

		for (size_t i = 0; i < nc; i++)
			rfs += rf(t,i)*ls[i];

		n[0] = gdt * -rfs.real();
		n[1] = gdt *  rfs.imag();
		n[2] = gdt * (gr(X,t) * lr[X] + gr(Y,t) * lr[Y] + gr(Z,t) * lr[Z] + b0m[pos] * TWOPI);

		RotationMatrix (n, rot);
		
		tmp[X] = rot[0]*ml[X] + rot[3]*ml[Y] + rot[6]*ml[Z];
		tmp[Y] = rot[1]*ml[X] + rot[4]*ml[Y] + rot[7]*ml[Z];
		tmp[Z] = rot[2]*ml[X] + rot[5]*ml[Y] + rot[8]*ml[Z];

		ml[X]  = tmp[0];
		ml[Y]  = tmp[1];
		ml[Z]  = tmp[2];

	}

	m(X,pos) = ml[X]; 
	m(Y,pos) = ml[Y]; 
	m(Z,pos) = ml[Z]; 

}


/**
 * @brief       Simulate single shot reception of freely precessing isochromat along gradient trajectory<br/>
 *              (i.e. forward Fourier transform incl. effect of Receive and b0 maps)<br/>
 *              Expects res to carry the correct size and dimensions
 *
 * @param  rxm  Receive sensitivities
 * @param  gr   Gradient trajectory
 * @param  r    Positions
 * @param  m0   Initial magnetisation state
 * @param  b0m  B0 map
 * @param  v    Verbose
 * @param  pos  Position of isochromat in r
 * @param  sig  Resulting signal (output) 
 */
void SimulateRecv (const Matrix<cplx>&   rxm, const Matrix<double>& gr, const Matrix<double>& r, const Matrix<double>& m0, 
				   const Matrix<double>& b0m, const double&         dt, const bool&           v, const size_t&        pos, 
				   const int&            tid,       Matrix<cplx>&  res) { 

	using namespace std;

	size_t nt = gr.Dim(1);      // Time points
	size_t nc = rxm.Dim(1);     // # channels

	Matrix<double> n   (3,1);   // Rotation axis
	Matrix<double> rot (3,3);   // Rotation matrix
	Matrix<double> m   (3,1);   // Magnetisation
	Matrix<double> tmp (3,1);   // Temporary magnetisation

	Matrix<double> lr  (3,1);   // Local spatial vector
	Matrix<cplx>   ls  (nc,1);  // Local sensitivity

	for (size_t c = 0; c < nc; c++) ls[c] = conj(rxm (pos,c));
	for (size_t i = 0; i <  3; i++) lr[i] =      r   (i,pos) ;

	// Starting magnetisation
	m[0] = m0(X,pos), m[1] = m0(Y,pos); m[2] = m0(Z,pos);

	double gdt = GAMMARAD * dt;

	// Run over time points
	for (size_t t = 0; t < nt; t++) {

		// Rotate magnetisation (only gradients)
		n[2] = gdt * (gr(X,t)*lr[X] + gr(Y,t)*lr[Y] + gr(Z,t)*lr[Z] + b0m[pos]*TWOPI);
		
		RotationMatrix (n, rot);
		
		tmp[0] = rot[0]*m[X] + rot[3]*m[Y] + rot[6]*m[Z];
		tmp[1] = rot[1]*m[X] + rot[4]*m[Y] + rot[7]*m[Z];
		tmp[2] = rot[2]*m[X] + rot[5]*m[Y] + rot[8]*m[Z];
		
		m[0]   = tmp[0];
		m[1]   = tmp[1];
		m[2]   = tmp[2];
		
		// Weighted contribution to all coils
		for (size_t c = 0; c < nc; c++)
			res(t,c,tid) += cplx(ls[c].real()*m[X], ls[c].imag()*m[Y]);
		
	}

}



/**
 * @brief Piece-wise constant bloch simulation
 *
 * @param  txm  Transmit sensitivity    (Nr  x Ntxc)
 * @param  rxm  Receive sensitivity     (Nr  x Nrxc) 
 * @param  rf   RF field                (Nt  x Nt  )
 * @param  gr   Gradient                (1-3 x Nt  )
 * @param  r    Spatial positions       (1-3 x Nr  )
 * @param  m0   Starting magnetisation  (3   x Nr  )
 * @param  b0m  B0 map                  (Nr        )
 * @param  v    Verbose                 (Scalar: false = only end, true = all time points)
 * @param  np   # parallel processes    (scalar)
 * 
 * @param  res  Result of simulation    (Nr  x 3 (x Nt))
 * @param  m    Resulting magnetisation (3   x Nt)
 */
void Simulate (const Matrix<cplx>&   txm, const Matrix<cplx>&   rxm, const Matrix<cplx>&   rf, 
			   const Matrix<double>&  gr, const Matrix<double>&   r, const Matrix<double>& m0, 
			   const Matrix<double>& b0m, const double&          dt, const bool&          exc, 
			   const bool&             v, const size_t&          np,
			   Matrix<cplx>&         res,       Matrix<double>&   m) {

	ticks  tic  = getticks();
	
	size_t nr   =   r.Dim(1);
	size_t ntxc = txm.Dim(1);
	size_t nrxc = rxm.Dim(1);
	size_t nt   =  gr.Dim(1);

	printf ("  Simulaing: %s on %04i isochromats ... ", (exc) ? "         excitation" : " signal acquisition", (int)nr); fflush(stdout);

	// Anything to do? ----------------

	if (gr.Size() < 1 || r.Size() < 1) {
		std::cout << "  Bailing out: %i Gradient step for %i isochromats? I don't think so!\n" << std::endl;
		return;
	}
	// --------------------------------

	Matrix<cplx> mres;

	if (!exc)
		mres = Matrix<cplx> (res.Dim(0), res.Dim(1), np);
	
#pragma omp parallel default (shared) 
	{
		
		omp_set_num_threads((int)np);
		int tid = omp_get_thread_num();
		
		if (exc) {
#pragma omp for schedule (guided, 10)
			for (size_t i = 0; i < nr; i++)
				SimulateExc (txm, rf, gr,  r, b0m, dt, v, i, tid, m);
		} else {
#pragma omp for  schedule (guided, 10)
			for (size_t i = 0; i < nr; i++)
				SimulateRecv (rxm, gr, r, m0, b0m, dt, v, i, tid, mres);
#pragma omp for  schedule (guided, 10)
			for (size_t i = 0; i < res.Size(); i++) {
				for (int p = 0; p < (int)np; p++)
					res[i] += mres[p*res.Size()+i];
				res[i] /= (float)nr;
			}
		}
 	}
	
	printf (" done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

}


