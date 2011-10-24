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
	
	double ar, ai, br, bi, hp, cp, sp, aiai, arai2, brbr,
		bibi, brbi2, arbi2, aibr2, arbr2, aibi2, phi;

	phi = sqrt(n[X]*n[X] + n[Y]*n[Y] + n[Z]*n[Z]);

	// Identity
	if (!phi)
		r = Matrix<double>::Id(3);
	
	else {
		// Cayley-Klein parameters
		hp =  phi/2;		
		cp =  cos(hp);
		sp =  sin(hp)/phi; /* /phi because n is unit length in defs. */
		ar =  cp;
		ai = -n[Z]*sp;
		br =  n[Y]*sp;
		bi = -n[X]*sp;
	
		 	/* Make auxiliary variables to speed this up	*/

		double arar  =   ar*ar;
		double aiai  =   ai*ai;
		double arai2 = 2*ar*ai;
		brbr  =   br*br;
		bibi  =   bi*bi;
		brbi2 = 2*br*bi;
		arbi2 = 2*ar*bi;
		aibr2 = 2*ai*br;
		arbr2 = 2*ar*br;
		aibi2 = 2*ai*bi;
		
		size_t i = 0;
		r[i++] =  arar  - aiai - brbr + bibi;
		r[i++] = -arai2 - brbi2;
		r[i++] = -arbr2 + aibi2;
		r[i++] =  arai2 - brbi2; 
		r[i++] =  arar  - aiai + brbr - bibi;
		r[i++] = -aibr2 - arbi2;
		r[i++] =  arbr2 + aibi2;
		r[i++] =  arbi2 - aibr2;
		r[i++] =  arar  + aiai - brbr - bibi;
		
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

	
	assert (txm.HDim() <= 3);
	assert (gr.Dim(1) == rf.Dim(1));
	
	size_t nt = gr.Dim(1);           // Time points
	size_t nc = txm.Dim(txm.HDim()); // # channels
	
	Matrix<double> n   (3,1);          // Rotation axis
	Matrix<double> rot (3,3);          // Rotation matrix
	Matrix<cplx>   ls  (nc,1);         // Local sensitivity
	Matrix<double> ml  (3,1);
	Matrix<double> tmp (3,1);

	ml[0] = 0.0; m[1] = 0.0; m[2] = 1.0;

	for (size_t i = 0; i < nc; i++)
		ls[i] = txm(pos,i);

	// Run over time points
	for (size_t t = 0; t < nt; t++) {

		cplx rfs = cplx (0.0,0.0);

		for (size_t i = 0; i < nc; i++)
			rfs += rf(t,i)*ls(pos,i);

		n[0] = -rfs.real();
		n[1] =  rfs.imag();
		n[2] = (gr(X,t) * r(X,pos) + gr(Y,t) * r(Y,pos) + gr(Z,t) * r(Z,pos) + b0m[pos] * TWOPI);

		n *= (GAMMARAD * dt);
		
		RotationMatrix (n, rot);
		
		tmp[0] = rot[0]*m[X] + rot[3]*m[Y] + rot[6]*m[Z];
		tmp[1] = rot[1]*m[X] + rot[4]*m[Y] + rot[7]*m[Z];
		tmp[2] = rot[2]*m[X] + rot[5]*m[Y] + rot[8]*m[Z];

		m[0] = tmp[0];
		m[1] = tmp[1];
		m[2] = tmp[2];
		
	}

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

	size_t nt = gr.Dim(1);           // Time points
	size_t nc = rxm.Dim(rxm.Dim(1)); // # channels

	Matrix<double> n   (3,1);          // Rotation axis
	Matrix<double> rot (3,3);          // Rotation matrix
	Matrix<double> m   (3,1);          // Magnetisation
	Matrix<double> tmp (3,1);          // Magnetisation
	Matrix<double> lr  (3,1);          // Local spatial vector

	Matrix<cplx>   ls  (nc,1);         // Local sensitivity

	for (size_t c = 0; c < nc; c++)
		ls[c] = conj(rxm (pos,c));

	for (size_t i = 0; i <  3; i++)
		lr[i] = r   (i,pos);

	// Starting magnetisation
	m[0] = m0(X,pos), m[1] = m0(Y,pos); m[2] = m0(Z,pos);

	// Run over time points
	for (size_t t = 0; t < nt; t++) {

		//printf ("%i\n", t);

		n[0] = 0.0;
		n[1] = 0.0;
		n[2] = GAMMARAD * dt * (gr(X,t) * lr[X] + gr(Y,t) * lr[Y] + gr(Z,t) * lr[Z] + b0m[pos] * TWOPI);

		RotationMatrix (n, rot);

		tmp[0] = rot[0]*m[X] + rot[3]*m[Y] + rot[6]*m[Z];
		tmp[1] = rot[1]*m[X] + rot[4]*m[Y] + rot[7]*m[Z];
		tmp[2] = rot[2]*m[X] + rot[5]*m[Y] + rot[8]*m[Z];

		m[0] = tmp[0];
		m[1] = tmp[1];
		m[2] = tmp[2];
		
		for (size_t c = 0; c < nc; c++)
			res(t,c,tid) += cplx(ls[c].real() * m[X], ls[c].imag() * m[Y]);

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
			   const Matrix<double>& b0m, const double&          dt, const bool           exc, 
			   const bool&             v, const size_t&          np,
			   Matrix<cplx>&         res,       Matrix<double>&   m) {

	ticks  tic  = getticks();
	
	size_t nr   =   r.Dim(1);
	size_t ntxc = txm.Dim(1);
	size_t nrxc = rxm.Dim(1);
	size_t nt   =  gr.Dim(1);

	printf ("  Simulaing Bloch equations for %i isochromats ... \n", (int)nr); fflush(stdout);

	// Anything to do? ----------------

	if (gr.Size() < 1 || r.Size() < 1) {
		std::cout << "  Bailing out: %i Gradient step for %i isochromats? I don't think so!\n" << std::endl;
		return;
	}
	// --------------------------------

	// Excitation? --------------------

	if (exc)
		assert(rf.Dim(1) == gr.Dim(1));

	printf ("  mode: %s\n", (exc) ? "excitation" : "acquisition");
	// --------------------------------
	
	
	Matrix<cplx> mres = Matrix<cplx>(res.Dim(0), res.Dim(1), np);
	
#pragma omp parallel default (shared) 
	{
		
		omp_set_num_threads((int)np);
		int tid = omp_get_thread_num();
		
		if (exc) {
			/*
			  #pragma omp for
			  for (size_t i = 0; i < nr; i++) 
			  SimulateExc (txm, rf, gr,  r, b0m, dt, v, i, tid, mres);
			*/
		} else {
#pragma omp for
			for (size_t i = 0; i < nr; i++)
				SimulateRecv (rxm, gr, r, m0, b0m, dt, v, i, tid, mres);
#pragma omp for schedule (dynamic, res.Size()/(int)np)
			for (size_t i = 0; i < res.Size(); i++) {
				for (int p = 0; p < (int)np; p++)
					res[i] += mres[p*res.Size()+i];
				res[i] /= (float)nr;
			}
		}
 	}
	
	printf ("  ... done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

}


