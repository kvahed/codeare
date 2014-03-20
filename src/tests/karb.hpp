#include "../modules/GradientTiming.hpp"
#include "../modules/SimulatedAnnealing.hpp"

bool karbtest (Connector& rc) {

	GradientParams gp;

	std::string in_file;// = std::string (base_dir + "/" + std::string(rc.GetElement("/config/data-in")->Attribute("fname")));
	std::string out_file;// = std::string (base_dir + "/" + rc.GetElement("/config/data-out")->Attribute("fname"));

	bool   simann  = false;
	size_t hpoints = 1;
	rc.Attribute ("simann",   &simann);
	rc.Attribute ("hpoints",  &hpoints);

	if (simann) {

		double coolrate, startt, finalt;
		size_t coolit;
		bool   verbose, exchange;
		
		rc.Attribute ("coolrate", &coolrate);
		rc.Attribute ("startt",   &startt);
		rc.Attribute ("finalt",   &finalt);
		rc.Attribute ("coolit",   &coolit);
		rc.Attribute ("verbose",  &verbose);
		rc.Attribute ("exchange", &exchange);
		
		SimulatedAnnealing sa (gp.k, coolit, startt, finalt, coolrate, verbose, exchange);
		sa.Cool();
		gp.k = sa.GetSolution();
		
	}
	
	rc.Attribute ("maxgrad", &(gp.mgr));
	rc.Attribute ("maxslew", &(gp.msr));
	rc.Attribute ("dt",      &(gp.dt));
	
	rc.Attribute ("gunits",  &(gp.gunits));
	rc.Attribute ("lunits",  &(gp.lunits));
	
	Matrix<double> x  = linspace<double> (0.0,1.0,size(gp.k,0));
	Matrix<double> xi = linspace<double> (0.0,1.0,size(gp.k,0)*hpoints);

	gp.k = interp1 (x, gp.k, xi, INTERP::AKIMA);

	printf ("\nComputing trajectory for ... \n");
	printf ("    [maxgrad: %.2f, maxslew: %.2f, dt: %.2e]\n\n", gp.mgr, gp.msr, gp.dt);
	
    SimpleTimer st ("VD spiral design");
	Solution s = ComputeGradient (gp);
    st.Stop();
	
    IOContext f = fopen (out_file.c_str(), WRITE);
    s.dump (f);
    fclose (f);
	

	return true;

}
