#include "../modules/GradientTiming.hpp"

template <ConType CT> bool slicegrad (Connector<CT>& rc) {

	GradientParams gp;
    
    IOContext fin (rc.GetElement("/config/data-in"), base, READ);
    gp.k = fin.Read<double>(rc.GetElement("/config/data-in/kin"));
    fin.Close();

	rc.Attribute ("maxgrad", &(gp.mgr));
	rc.Attribute ("maxslew", &(gp.msr));
	rc.Attribute ("dt",      &(gp.dt));
	
	rc.Attribute ("gunits",  &(gp.gunits));
	rc.Attribute ("lunits",  &(gp.lunits));

	//Matrix<double> x  = linspace<double> (0.0,1.0,size(gp.k,0));
	//Matrix<double> xi = linspace<double> (0.0,1.0,size(gp.k,0)*hpoints);
	//gp.k = k;//interp1 (x, gp.k, xi, INTERP::AKIMA);

    Solution s = ComputeGradient (gp);

    IOContext fout (rc.GetElement("/config/data-out"), base, WRITE);
    s.dump (fout);
    fout.Close();

    return true;
}
