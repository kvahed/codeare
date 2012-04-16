#include "../modules/GradientTiming.hpp"

template <class T> bool
karbtest (Connector<T>* rc) {

	GradientParams gp;

	std::string    cf  = std::string (base + std::string(config));
	rc->ReadConfig (cf.c_str());

#ifdef HAVE_MAT_H

	// Gradients
	Read    (gp.k, rc->GetElement("/config/data/k"),    base);

#endif

	gp.k /= 10.0;

	size_t nn = size(gp.k,0);
	double nd = (double) nn;

	/*	Matrix<double> xx = linspace<double> (0, nd, nn);
	Matrix<double> XX = linspace<double> (0, nd, nn*100);

	MXDump (xx, "x.mat", "x");
	MXDump (XX, "xi.mat", "xi");


	gp.k = interp1 (xx, gp.k, XX);*/
	
	MXDump (gp.k, "gpk.mat", "gpk");

	rc->Attribute ("maxgrad", &(gp.mgr));
	rc->Attribute ("maxslew", &(gp.msr));
	rc->Attribute ("dt",      &(gp.dt));

	rc->Attribute ("gunits",  &(gp.gunits));
	rc->Attribute ("lunits",  &(gp.lunits));
	
	printf ("Computing trajectory for ... \n");
	printf ("    [maxgrad: %.2f, maxslew: %.2f]\n", gp.mgr, gp.msr);
	
	ticks start = getticks();
	Solution s = ComputeGradient (gp);
	printf ("... done. WTime: %.4f seconds.\n\n", elapsed(getticks(), start) / Toolbox::Instance()->ClockRate());
	
#ifdef HAVE_MAT_H	
	
	MATFile* mf = matOpen ("spiral.mat", "w");
	
	if (mf == NULL) {
		printf ("Error creating file %s\n", "");
		return false;
	}
	
	MXDump (s.g, mf, "g");
	MXDump (s.k, mf, "k");
	MXDump (s.s, mf, "sr");
	MXDump (s.t, mf, "t");
	
	if (matClose(mf) != 0) {
		printf ("Error closing file %s\n", "");
		return false;
	}
	
#endif
	
	return true;

}
