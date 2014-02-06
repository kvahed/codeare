
/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
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

	double shots;       /**< @brief # shots */
	double res;         /**< @brief Maximum resolution */
	Matrix<double> fov; /**< @brief FOV vector */
	Matrix<double> rad; /**< @brief Corresponding radius vector */
	double mgr;         /**< @brief G max */
	double msr;         /**< @brief Slew max */
	double dt;          /**< @brief Sampling duration (i.e. delta t) */
	int    gunits;
	int    lunits;     

};


Solution VDSpiral (SpiralParams& sp) {

	GradientParams gp;

	Matrix<double>& fov = sp.fov; 
	Matrix<double>& rad = sp.rad; 
	double k_max, fov_max, dr;
	Matrix<double> r, theta;
	long n = 0;

	assert (numel(rad) >= 2);
	assert (isvec(rad) == isvec(fov));
	assert (numel(rad) == numel(fov));

	k_max   = 5.0 / sp.res;
	fov_max = m_max(fov);

	dr  = sp.shots / (fov_max);
	n   = size(fov,1)*100;
	r   = linspace<double> (0.0, k_max, n);
	
	Matrix<double> x = k_max*rad;
	fov = interp1 (x, fov, r, INTERP::LINEAR);

	dr  = sp.shots / (1500.0 * fov_max);
	n   = ceil (k_max/dr);
	x   = r;
	r   = linspace<double> (0.0, k_max, n);

	fov = interp1 (x, fov, r, INTERP::AKIMA);

	theta = cumsum ((2 * PI * dr / sp.shots) * fov);

	gp.k = Matrix<double> (numel(r), 3);

	for (size_t i = 0; i < numel(r); i++) {
		gp.k(i,0) = r[i] * cos (theta[i]);
		gp.k(i,1) = r[i] * sin (theta[i]);
	}

	gp.mgr     = sp.mgr;
	gp.msr     = sp.msr;
	gp.dt      = sp.dt;
	gp.gunits  = sp.gunits;
	gp.lunits  = sp.lunits;

	Solution s = ComputeGradient (gp);
	
	return s;

}

