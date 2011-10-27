#include "BlochSim.hpp"

/**
 * @brief         Time reversal for RF 
 *
 * @param  signal Acquired signal
 * @param  jac    Jacobian determinant j(k) i.e. density compensation
 * @param  pulse  Excitation pulse(s)
 */
void TimeReverseRF (const Matrix<cplx>& signal, const Matrix<double>& jac, Matrix<cplx>& pulse) {
	
	size_t nt = signal.Dim(0);
	size_t nc = signal.Dim(1);
	
	for (size_t i = 0; i < nt; i++)
		for (size_t c = 0; c < nc; c++)
			pulse(i,c) = (signal(nt-1-i,c)*(float)jac[nt-1-i]);
	
}


/**
 * @brief         Time reversal for gradient trajectory
 * 
 * @param  acqgr  Acquisition gradients
 * @param  excgr  Excitation gradients
 */
void TimeReverseGR (const Matrix<double>& acqgr, Matrix<double>& excgr) {
	
	size_t nt = acqgr.Dim(1);
	
	for (size_t i = 0; i < nt; i++) {
		
		excgr(0,i) = -acqgr(0,nt-1-i); 
		excgr(1,i) = -acqgr(1,nt-1-i); 
		excgr(2,i) = -acqgr(2,nt-1-i);
		
	}
	
}


/**
 * @brief Simulation mode
 */
enum sim_mode {
	ACQUIRE, 
	EXCITE
};

/**
 * @brief          CG NR algorithm
 * 
 * @param  maxit   Maximum # CG iterations
 * @param  eps     Convergence criterium 
 * @param  rxm     Receive sensitivities (Limited to pattern)
 * @param  txm     Transmit sensitivities (Limited to "care-region") in most cases field od view
 * @param  acqg    Acquisition gradients
 * @param  rv      Spatial vectors for target
 * @param  sv      Spatial vectors for sample
 * @param  target  Target magnetisation
 * @param  sample  Excitation sample (i.e. Care region)
 * @param  rb0     Receive b0
 * @param  sb0     Transmit b0 (in general, both b0 identical with different signs yet 
 *                              each limited to appropriate size)
 * @param  dt      Delta t. Simualtion time interval
 * @param  verb    Verbosity
 * @param  np      # processing threads
 * @param  mag     Resulting magnetisations from excitation simulations (Outgoing)
 * @param  rf      Acquired RF pulse shapes (Outgoing)
 */
/*
void CGNR (//IN
		   const int&          maxit, const float&             eps, const Matrix<cplx>&   rxm, 
		   const Matrix<cplx>&   txm, const Matrix<double>&   acqg, const Matrix<double>&  rv, 
		   const Matrix<double>&  sv, const Matrix<double>& target, const Matrix<double>& rb0, 
		   const Matrix<double>& sb0, const Matrix<double>&    jac, const float&           dt, 
		   const bool&          verb, const int&                np,
		   // OUT
		   Matrix<cplx>          mag, Matrix<cplx>              rf                               
		   ) {
	
	int             iters = 0;
	cplx               rn = 0.0;
	float              an = 0.0;
	cplx             rtmp = cplx(0.0,0.0);
	
	vector<float> residue ;
	Matrix<double>   excg   (acqg.Dim(0),acqg.Dim(1));
	Matrix<cplx>      res   (acqg.Dim(1), txm.Dim(1));
	
	TimeReverseGR (acqg, excg);
	
	Matrix<cplx> a(m_N[0], m_N[1], m_N[2]);
	Matrix<cplx> s(m_M, m_Nc);

	// Convergence loop
	for (iters = 0; iters < maxit; iters++) {
		
		rn = pow(r.Norm().real(), 2);
		residue.push_back(rn/an);
		printf ("  %03i: CG residuum: %.9f\n", iters, residue.at(iters));
		
		// Convergence ? ----------------------------------------------
		if (std::isnan(residue.at(iters)) || residue.at(iters) <= eps)
			break;
		
		// E^H*E
		Simulate (txm, rxm, rf, acqg, tv, target, rb0, dt, ACQUIRE, verb, np, res, mag); // Simulate Bloch receive mode
		TimeReverseRF (res, jac, rf);                                                    // Time reversal
		Simulate (txp, rxm, rf, excg, sv, sample, sb0, dt,  EXCITE, verb, np, res, mag); // Simulate Bloch transmit mode
		
		rtmp    = (rn / (target.dotc(mag)));
		a      += (target    * rtmp);
		s      += (stmp * rtmp);
		r      -= (mag    * rtmp);
		target *= cplx(pow(r.Norm().real(), 2)/rn);
		target += r;
		
		// Verbose out put keeps all intermediate steps ---------------
		if (m_verbose) {
			image.Expand(m_dim);  
			signals.Expand(2);
			memcpy (  &image[(iters+1)*a.Size()], &a[0], a.Size()*sizeof(cplx));
			memcpy (&signals[(iters+1)*s.Size()], &s[0], s.Size()*sizeof(cplx));
		}
		
	}	
	
}
*/
