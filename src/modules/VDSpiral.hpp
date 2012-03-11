
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

#include "GradientTiming.hpp"

/**
 * @brief 2D Spiral parameters
 */
struct SpiralParams {

	size_t shots;       /**< @brief # shots */
	double res;         /**< @brief Maximum resolution */
	Matrix<double> fov; /**< @brief FOV vector */
	Matrix<double> rad; /**< @brief Corresponding radius vector */
	double mgr;         /**< @brief G max */
	double msr;         /**< @brief Slew max */
	double dt;          /**< @brief Sampling duration (i.e. delta t) */

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


Spiral VDSpiral (SpiralParams& sp) {

	Spiral spir;
	GradientParams gp;

	Matrix<double>& fov = sp.fov; 
	Matrix<double>& rad = sp.rad; 
	double k_max, fov_max, dr;
	Matrix<double> r, theta;
	long n = 0;

	assert (rad.Size() >= 2);
	assert (Is1D(rad) == Is1D(fov));
	assert (rad.Size() == fov.Size());

	k_max   = 5.0 / sp.res;
	fov_max = max(fov);

	dr  = ((double) sp.shots) / (fov_max);
	n   = size(fov,1)*100;
	r   = Matrix<double>::LinSpace(0.0, k_max, n);
	
	Matrix<double> x = k_max*rad;
	fov = interp1 (x, fov, r, INTERP::LINEAR);

	dr  = ((double) sp.shots) / (1500.0 * fov_max);
	n   = ceil (k_max/dr);
	x   = r;
	r   = Matrix<double>::LinSpace(0.0, k_max, n);

	fov = interp1 (x, fov, r, INTERP::AKIMA);

	MXDump (fov, "fov.mat", "fov");

	theta = cumsum((2 * PI * dr / sp.shots) * fov);

	gp.k = Matrix<double> (numel(r),3);

	for (size_t i = 0; i < numel(r); i++) {
		gp.k(i,0) = r[i] * cos (theta[i]);
		gp.k(i,1) = r[i] * sin (theta[i]);
	}

	gp.mgr = sp.mgr;
	gp.msr = sp.msr;
	gp.dt  = sp.dt;
	
	Solution s = ComputeGradient (gp);
	
	spir.k = s.k;
	spir.g = s.g;
	spir.s = s.s;
	spir.t = s.t;

	return spir;

}

