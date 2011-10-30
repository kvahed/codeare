#include "CPUSimulator.hpp"

using namespace RRStrategy;

#define GAMMARAD 2.6753e3
#define TWOPI	 6.283185


/**
 * @brief       Rotation matrix around vector n
 *
 * @param  n    Rotation axis
 * @param  r    Rotation matrix (output)
 */
void 
Rotate (const Matrix<double>& n, Matrix<double>& r) {
	
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


CPUSimulator::CPUSimulator (SimulationBundle& sb) {
	
	m_sb  = sb;
	m_sig = NEW (Matrix<cplx> (sb.rf->Dim(0), sb.rf->Dim(0), sb.np));

}



void 
CPUSimulator::SimulateExc  (const size_t& pos) {

	
	size_t nt = m_sb.egr->Dim(1);
	assert (nt == m_sb.rf->Dim(0));
	size_t nc = m_sb.sb1->Dim(1); // # channels
	
	Matrix<double> n   ( 3,1);  // Rotation axis
	Matrix<double> rot ( 3,3);  // Rotation matrix
	Matrix<double> ml  ( 3,1);  // Local magnetisation
	Matrix<double> tmp ( 3,1);  // Local magnetisation

	Matrix<double> lr  (3,1);    // Local spatial vector
	Matrix<cplx>   ls  (nc,1);   // Local sensitivity

	for (size_t i = 0; i <  3; i++)
		lr[i] = m_sb.sr->At(i,pos);
	for (size_t i = 0; i < nc; i++)
		ls[i] = m_sb.sb1->At(pos,i);

	ml[Z] = 1.0;

	double gdt = GAMMARAD * m_sb.dt;
	
	// Run over time points
	for (size_t t = 0; t < nt; t++) {
		
		size_t rt = nt-1-t;

		cplx rfs = cplx (0.0,0.0);

		for (size_t i = 0; i < nc; i++)
			rfs += m_sb.rf->At(rt,i)*ls[i];

		rfs *= m_sb.jac->At(rt);

		n[0] = gdt * -rfs.real();
		n[1] = gdt *  rfs.imag();
		n[2] = gdt * (m_sb.egr->At(X,rt) * lr[X] + m_sb.egr->At(Y,rt) * lr[Y] + m_sb.egr->At(Z,rt) * lr[Z] + m_sb.sb0->At(pos) * TWOPI);

		Rotate (n, rot);
		
		tmp[X] = rot[0]*ml[X] + rot[3]*ml[Y] + rot[6]*ml[Z];
		tmp[Y] = rot[1]*ml[X] + rot[4]*ml[Y] + rot[7]*ml[Z];
		tmp[Z] = rot[2]*ml[X] + rot[5]*ml[Y] + rot[8]*ml[Z];

		ml[X]  = tmp[0];
		ml[Y]  = tmp[1];
		ml[Z]  = tmp[2];

	}

	m_sb.magn->At(X,pos) = ml[X];
	m_sb.magn->At(Y,pos) = ml[Y];
	m_sb.magn->At(Z,pos) = ml[Z]; 

}


void 
CPUSimulator::SimulateRecv (const size_t& pos, const int& tid) { 

	using namespace std;

	size_t nt = m_sb.agr->Dim(1);      // Time points
	size_t nc = m_sb.tb1->Dim(1);     // # channels

	Matrix<double> n   (3,1);   // Rotation axis
	Matrix<double> rot (3,3);   // Rotation matrix
	Matrix<double> m   (3,1);   // Magnetisation
	Matrix<double> tmp (3,1);   // Temporary magnetisation

	Matrix<double> lr  (3,1);   // Local spatial vector
	Matrix<cplx>   ls  (nc,1);  // Local sensitivity

	for (size_t c = 0; c < nc; c++) ls[c] = conj(m_sb.tb1->At(pos,c));
	for (size_t i = 0; i <  3; i++) lr[i] =      m_sb.tr->At(i,pos) ;

	// Starting magnetisation
	m[0] = m_sb.tm->At(X,pos), m[1] = m_sb.tm->At(Y,pos); m[2] = m_sb.tm->At(Z,pos);

	double gdt = GAMMARAD * m_sb.dt;

	// Run over time points
	for (size_t t = 0; t < nt; t++) {

		// Rotate magnetisation (only gradients)
		n[2] = gdt * (m_sb.agr->At(X,t)*lr[X] + m_sb.agr->At(Y,t)*lr[Y] + m_sb.agr->At(Z,t)*lr[Z] + m_sb.tb0->At(pos)*TWOPI);
		
		Rotate (n, rot);
		
		tmp[0] = rot[0]*m[X] + rot[3]*m[Y] + rot[6]*m[Z];
		tmp[1] = rot[1]*m[X] + rot[4]*m[Y] + rot[7]*m[Z];
		tmp[2] = rot[2]*m[X] + rot[5]*m[Y] + rot[8]*m[Z];
		
		m[0]   = tmp[0];
		m[1]   = tmp[1];
		m[2]   = tmp[2];
		
		// Weighted contribution to all coils
		for (size_t c = 0; c < nc; c++)
			m_sig->At(t,c,tid) += cplx(ls[c].real()*m[X], ls[c].imag()*m[Y]);
		
	}

}



void 
CPUSimulator::Simulate () {

	this->Simulate (ACQUIRE);

	m_sb.Dump("sb.mat");

	for (size_t i = 0; i < m_sb.rf->Dim(0); i++)
		for (size_t j = 0; j < m_sb.rf->Dim(1); j++)
			m_sb.rf->At(i,j) *= m_sb.jac->At(i);


	//this->Simulate (EXCITE);

}


void
CPUSimulator::Simulate (const bool& mode) {
	
	ticks            tic  = getticks();       // Start timing
	size_t           nr   = (mode) ? m_sb.sr->Dim(1) : m_sb.tr->Dim(1);  // Spatial vectors (target)
	size_t           ntxc = m_sb.sb1->Dim(1); // transmit channels
	size_t           nrxc = m_sb.tb1->Dim(1); // receive channels
	size_t           nt   = m_sb.agr->Dim(1); // time points
	size_t           ns   = m_sb.rf->Size();  // signals samples
	
	// Expect same #ch
	assert (ntxc == nrxc);

	printf ("  Simulating: %02ich %s on %04i isochromats ... ", (int)ntxc, (mode) ? "         excitation" : " signal acquisition", (int)nr); fflush(stdout);

	// Anything to do? ----------------

	if (nt < 1 || nr < 1) {
		printf ("    %i Gradient step for %i isochromats? I don't think so!\n", nt, nr);
		return;
	}
	// --------------------------------

#pragma omp parallel default (shared) 
	{
		
		omp_set_num_threads(m_sb.np);
		int tid = omp_get_thread_num();
		
		if (mode) {
#pragma omp for schedule (guided, 10)
			for (size_t i = 0; i < nr; i++)
				SimulateExc (i);
		} else {
#pragma omp for  schedule (guided, 10)
			for (size_t i = 0; i < nr; i++)
				SimulateRecv (i, tid);
#pragma omp for  schedule (guided, 10)
			for (size_t i = 0; i < ns; i++) {
				m_sb.rf->At(i) = cplx(0.0,0.0);
				for (int p = 1; p < m_sb.np; p++)
					m_sb.rf->At(i) += m_sig->At(p*ns+i);
				m_sb.rf->At(i) /= (float)nr;
			}
		}
	}

	printf (" done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());

}

