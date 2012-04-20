#include "../modules/GradientTiming.hpp"
#include "../modules/SimulatedAnnealing.hpp"

template <class T> bool
karbtest (Connector<T>* rc) {

	GradientParams gp;

	std::string    cf  = std::string (base + std::string(config));
	rc->ReadConfig (cf.c_str());

#ifdef HAVE_MAT_H

	// Gradients
	Read    (gp.k, rc->GetElement("/config/data/k"),    base);

#endif

	bool   simann  = false;
	size_t hpoints = 1;
	rc->Attribute ("simann",   &simann);
	rc->Attribute ("hpoints",  &hpoints);

	if (simann) {

		double coolrate, startt, finalt;
		size_t coolit;
		bool   verbose, exchange;
		
		rc->Attribute ("coolrate", &coolrate);
		rc->Attribute ("startt",   &startt);
		rc->Attribute ("finalt",   &finalt);
		rc->Attribute ("coolit",   &coolit);
		rc->Attribute ("verbose",  &verbose);
		rc->Attribute ("exchange", &exchange);
		
		SimulatedAnnealing sa (gp.k, coolit, startt, finalt, coolrate, verbose, exchange);
		sa.Cool();
		gp.k = sa.GetSolution();
		
	}
	
	rc->Attribute ("maxgrad", &(gp.mgr));
	rc->Attribute ("maxslew", &(gp.msr));
	rc->Attribute ("dt",      &(gp.dt));
	
	rc->Attribute ("gunits",  &(gp.gunits));
	rc->Attribute ("lunits",  &(gp.lunits));
	
	
	Matrix<double> x  = linspace<double> (0.0,1.0,size(gp.k,0));
	Matrix<double> xi = linspace<double> (0.0,1.0,size(gp.k,0)*hpoints);

	gp.k = interp1 (x, gp.k, xi, INTERP::AKIMA);

	printf ("Computing trajectory for ... \n");
	printf ("    [maxgrad: %.2f, maxslew: %.2f, dt: %.2e]\n", gp.mgr, gp.msr, gp.dt);
	
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
