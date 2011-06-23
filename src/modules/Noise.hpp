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
AddPseudoRandomNoise (Matrix<raw>* m, float max) {

	SetSeed();

	for (int i = 0; i < m->Size(); i++)
		m->at(i) = m->at(i) + max*WhiteNoise();
	
	return OK;

}


