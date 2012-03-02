#include "Algos.hpp"
#include "CX.hpp"
#include "Math.hpp"
#include "Interpolate.hpp"


#define M_GAMMA 4.257

struct Solution {
	
	Matrix<double> k;
	Matrix<double> g;
	Matrix<double> s;
	double t;
	
};


struct GradientParams {
	
	Matrix<double> k;  // K-Space trajectory [1/cm]
	double         mgr;  // Maximum allowed gradient
	double         msr;  // Maximum allowed slewrate
	double         dt;   // Sampling interval
	
};


inline double RugeKutta (const double& s, const double& ds, const double& st, 
						 const Matrix<double>& k, const double& smax, const double& L) {

	double k1, k2, k3, k4;
	Matrix<double> ak = abs(k) ^ 2;

	double pgs = pow(M_GAMMA*smax,2);

	k1 = ds * sqrt (pgs - ak[0] * pow (st       ,4)) /  st;
	k2 = ds * sqrt (pgs - ak[1] * pow (st + k1/2,4)) / (st + k1/2);
	k3 = ds * sqrt (pgs - ak[1] * pow (st + k2/2,4)) / (st + k2/2);
	k4 = ds * sqrt (pgs - ak[2] * pow (st + k3/2,4)) / (st + k3/2);

	return k1/6 + k2/3 + k3/3 + k4/6;

}


struct SDIn {

	PolyVal *pkx, *pky, *pkz;
	Matrix<double> pos, s;
	double mgr, msr;

	~SDIn () {
		
		delete pkx;
		delete pky;
		delete pkz;
		
	}

};


struct SDOut {

	Matrix<double> k;
	Matrix<double> phi;
	
};



SDOut SDMax (SDIn& sdin) {

	SDOut so;
	size_t ssp;

	Matrix<double> dpp (ssp,1), dpm (ssp,1);

	for (size_t i = 0; i < ssp-1; i++)
		dpp[i] = sdin.pos[i+1] = sdin.pos[i]; 
	dpp[ssp-1] = dpp[ssp-2];

	MXDump (dpp, "dpp.mat", "dpp");

	

	return so;

}
/*
  function [sdot, k] = sdotMax(PP, p_of_s, s, gmax, smax)

  gamma = 4.257;
  
  s = s(:);
  dp_p = p_of_s([2:end,end]) - p_of_s; , dp_p(end) = dp_p(end-1);
  dp_m = p_of_s - p_of_s([1,1:end-1]);, dp_m(1) = dp_m(2);
  ds_p = s([2:end,end]) - s; , ds_p(end) = ds_p(end-1);
  ds_m = s - s([1,1:end-1]);, ds_m(1) = ds_m(2);
  
  Cs_p = (ppval(PP,p_of_s + dp_p) - ppval(PP, p_of_s))./ds_p;
  Cs_m = (ppval(PP,p_of_s) - ppval(PP, p_of_s-dp_m))./ds_m;
  Cs = Cs_p/2 + Cs_m/2;
  Css = (Cs_p - Cs_m)./(ds_m/2+ds_p/2);
  k = abs(Css);
  % fix edge numerical problems
  k(end) = k(end-1);
  k(1) = k(2);
  
  % calc I constraint curve (maximum gradient)
  sdot1 = gamma*gmax*ones(size(s));
  
  % calc II constraint curve (curve curvature dependent)
  sdot2 = sqrt(gamma*smax ./ (abs(k)+eps));
  
  % calc total constraint
  sdot = min([sdot1, sdot2],[],2);
*/



	


/**
 * @brief       Compute real world k-space / gradient solution with hardware limits
 *
 * @param  p    Parameters
 * @return      Solution
 */
Solution ComputeGradient (GradientParams& gp) {

	Solution s;
	size_t ups = 100, ts = 0;

	Matrix<double> op, np, sop;

	op = Matrix<double>::LinSpace (0.0, (double)(size(gp.k,0)-1), size(gp.k,0));
	np = Matrix<double>::LinSpace (0.0, (double)(size(gp.k,0)-1), (size(gp.k,0)-1)*ups + 1);

	s.k = interp1 (op, gp.k, np, INTERP::CSPLINE);
	size_t snp = size(s.k,0);
	s.g = Matrix<double>(snp,3);
	sop = Matrix<double>(snp,1);

	SDIn sdin;

	sdin.pkx = new PolyVal (np, (double*)&(s.k)[0],     INTERP::CSPLINE);
	sdin.pky = new PolyVal (np, (double*)&(s.k)[1*snp], INTERP::CSPLINE);
	sdin.pkz = new PolyVal (np, (double*)&(s.k)[2*snp], INTERP::CSPLINE);

	sdin.mgr = gp.mgr;
	sdin.msr = gp.msr;

	for (size_t i = 0; i < snp-1; i++) {
		s.g (i,0) = s.k(i+1,0) - s.k(i,0);
		s.g (i,1) = s.k(i+1,1) - s.k(i,1);
		s.g (i,2) = s.k(i+1,2) - s.k(i,2);
	}

	for (size_t i = 0; i < 3; i++)
		s.g (snp-1,i) = s.g (snp-2,i); 

	s.g *= (double)ups;

	sop [0] = 0.0;

	for (size_t i = 0; i < snp-1; i++)
		sop [i+1] = sop[i] + sqrt (pow(s.g(i,0),2) + pow(s.g(i,1),2) + pow(s.g(i,2),2));

	sop /= (double)ups;

	double st0 = M_GAMMA * gp.msr * gp.dt;
	double ds  = st0 * gp.dt / 3.0;
	
	double L = max(sop);
	ts = ceil(L/ds);

	Matrix<double> sh, sta;

	sdin.s = Matrix<double>::LinSpace (0.0, L,   ts);
	sh     = Matrix<double>::LinSpace (0.0, L, 2*ts);
	sta    = Matrix<double>::Zeros (size(sdin.s,0), 1);

	sta[0] = 0.0;

	Matrix<double> posh = interp1 (sop, np, sh, INTERP::CSPLINE);
	Matrix<double> pos  (size(posh,0)/2,1);

	for (size_t i = 0; i < size(pos,1); i++)
		pos [i] = posh[i*2];

	SDOut sdout = SDMax (sdin);

	return s;
	
}


/**
 * @brief 2D Spiral parameters
 */
struct SpiralParams {

	size_t shots;       // # shots
	double res;         // Maximum resolution
	Matrix<double>* fov; // FOV vector
	Matrix<double>* rad; // Corresponding radius vector
	double mgr;         /**< @brief G max */
	double msr;         /**< @brief Slew max */
	double dt;          /**< @brief Sampling duration (i.e. delta t)*/

};


/**
 * @brief 2D Spiral parameters
 */
struct Spiral {

	Matrix<double> k;    /**< @brief k-space trajectory */
	Matrix<double> g;    /**< @brief Gradient amplitudes */
	Matrix<double> s;    /**< @brief Slew rate */
	Matrix<double> t;    /**< @brief Time */

};


Spiral VDSpiral (const SpiralParams& sp) {

	Spiral spir;
	GradientParams gp;

	Matrix<double> &fov = *(sp.fov); 
	Matrix<double> &rad = *(sp.rad); 
	double k_max, fov_max, dr;
	Matrix<double> r, theta;
	long n = 0;

	assert (rad.Size() >= 2);
	assert (Is1D(rad) == Is1D(fov));
	assert (rad.Size() == fov.Size());

	k_max   = 5.0 / sp.res;
	fov_max = max(fov);
	printf ("  max(fov) %f\n", fov_max);

	dr  = ((double) sp.shots) / (10.0 * fov_max);
	n   = ceil (k_max/dr);
	r   = Matrix<double>::LinSpace(0.0, k_max, n);
	
	Matrix<double> x = k_max*rad;
	fov = interp1 (x, fov, r, INTERP::LINEAR);
	
	dr  = ((double) sp.shots) / (1500.0 * fov_max);
	n   = ceil (k_max/dr);
	x   = r;
	r   = Matrix<double>::LinSpace(0.0, k_max, n);

	fov = interp1 (x, fov, r, INTERP::AKIMA);
	
	theta = cumsum((2 * PI * dr / sp.shots) * fov);

	gp.k = Matrix<double> (numel(r),3);

	for (size_t i = 0; i < numel(r); i++) {
		gp.k(i,0) = r[i] * cos (theta[i]);
		gp.k(i,1) = r[i] * sin (theta[i]);
	}

	gp.mgr = sp.mgr;
	gp.msr = sp.msr;
	gp.dt  = sp.dt;
	
	MXDump (r,     "r.mat",   "r");
	MXDump (theta, "t.mat",   "t");

	Solution s = ComputeGradient (gp);

	return spir;

}

