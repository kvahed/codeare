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

#define M 714025
#define IA 1366
#define IC 150889
#define LM 2147483647
#define LAM (1.0/LM)
#define LA 16807
#define LR 2836
#define LQ 127773

long idum;

void 
SetSeed (double i=1349555.0) { 

	idum = (long)i;

}


double 
Uniform () {
	
	double  r1;
	long    hi;

	hi      = idum/LQ;
	idum    = LA*(idum-hi*LQ) - LR*hi;
	
	if (idum < 0) 
		idum += LM;
	
	r1  = LAM*idum;

	return(r1);

}


std::complex<float> 
WhiteNoise () {

	float fac, r, v1, v2;

	do {
		v1 = (2.0 * Uniform()) - 1.0;
		v2 = (2.0 * Uniform()) - 1.0;
		r = (v1*v1) + (v2*v2);
	} while (r >= 1.0);

	fac = sqrt(-2.0 * log(r) / r);

	return ( raw( v2*fac, v1*fac ) );

}


RRSModule::error_code
AddPseudoRandomNoise (Matrix<raw>& m, const float& max) {

	SetSeed();

	for (int i = 0; i < m.Size(); i++)
		m[i] = m[i] + max * WhiteNoise();
	
	return OK;

}


