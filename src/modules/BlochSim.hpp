#include "OMP.hpp"

#define GAMMARAD 26753.0
#define TWOPI	     6.283185

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
	
	double ar, ai, br, bi, hp, cp, sp, arar, aiai, arai2, brbr, 
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

		arar  =   ar*ar;
		aiai  =   ai*ai;
		arai2 = 2*ar*ai;
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
				   Matrix<cplx>&         res) { 

	assert (rxm.HDim() <= 3);

	size_t nt = gr.Dim(1);           // Time points
	size_t nc = rxm.Dim(rxm.HDim()); // # channels

	Matrix<double> n   (3,1);          // Rotation axis
	Matrix<double> rot (3,3);          // Rotation matrix
	Matrix<double> m   (3,1);          // Magnetisation
	Matrix<cplx>   ls  (nc,1);         // Local sensitivity
	Matrix<double> lr  (3,1);          // Local spatial vector

	for (size_t c = 0; c < nc; c++)
		ls(c) = rxm (pos,c);

	for (size_t i = 0; i <  3; i++)
		lr(i) = r   (i,pos);

	// Starting magnetisation
	m = m0;
	
	// Run over time points
	for (size_t t = 0; t < nt; t++) {

		n[0] = 0.0;
		n[1] = 0.0;
		n[2] = GAMMARAD * dt * (gr(X,t) * lr(X) + gr(Y,t) * lr(Y) + gr(Z,t) * lr(Z) + b0m(pos) * TWOPI);
		
		RotationMatrix (n, rot);
		
		m = rot->*m;

		for (size_t c = 0; c < nc; c++) {

			res(X,t,c) = ls(c).real() * m(X);
			res(Y,t,c) = ls(c).imag() * m(Y);
			res(Z,t,c) =                m(Z);

		}

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
void SimulateExc  (const Matrix<cplx>&   txm, const Matrix<cplx>& rf, const Matrix<double>& gr, const Matrix<double>& r, 
				   const Matrix<double>& b0m, const double&       dt, const bool&            v, const size_t&       pos, 
				   Matrix<double>& m) {

	
	assert (txm.HDim() <= 3);
	assert (gr.Dim(1) == rf.Dim(1));
	
	size_t nt = gr.Dim(1);           // Time points
	size_t nc = txm.Dim(txm.HDim()); // # channels
	
	Matrix<double> n   (3,1);          // Rotation axis
	Matrix<double> rot (3,3);          // Rotation matrix
	Matrix<cplx>   ls  (nc,1);         // Local sensitivity

	for (size_t i = 0; i < nc; i++)
		ls(i) = txm(pos,i);

	// Run over time points
	for (size_t t = 0; t < nt; t++) {

		cplx rfs = cplx (0.0,0.0);

		for (size_t i = 0; i < nc; i++)
			rfs += rf(t,i)*ls(pos,i);

		n[0] = - rfs.real();
		n[1] = + rfs.imag();
		n[2] =   (gr(X,t) * r(X,pos) + gr(Y,t) * r(Y,pos) + gr(Z,t) * r(Z,pos) + b0m(pos) * TWOPI);

		n *= (GAMMARAD * dt);
		
		RotationMatrix (n, rot);
		
		m = rot->*m;

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

	size_t nr   =   r.Dim(1);
	size_t ntxc = txm.Dim(1);
	size_t nrxc = rxm.Dim(1);
	size_t nt   =  gr.Dim(1);

	// Anything to do?
	if (gr.Size() < 1 || r.Size() < 1) return;

	// Excitation?
	if (exc)
		assert(rf.Dim(1) == gr.Dim(1));

#pragma omp parallel default (shared) 
	{
		
		omp_set_num_threads(np);
		int tid = omp_get_thread_num();
		
		if (exc) {

			res.Dim(0) = nt;
			res.Dim(1) = nrxc;
			res.Reset();
			
#pragma omp for
			for (size_t i = 0; i < nr; i++) 
				SimulateExc (txm, rf, gr, r, b0m, dt, v, i, m);
			
		} else {

			m.Dim(0) = 3;
			m.Dim(1) = nr;
			m.Reset();
			
#pragma omp for
			for (size_t i = 0; i < nr; i++) 
				SimulateRecv (rxm, gr, r, m0, b0m, dt, v, i, res);
			
		}
		
	}
	
}


