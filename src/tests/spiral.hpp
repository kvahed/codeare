#include "../modules/VDSpiral.hpp"
#include "Statistics.hpp"
#include "Tokenizer.hpp"

template <class T> bool
vdspiraltest (Connector<T>* rc) {

	SpiralParams p;

	std::string    cf  = std::string (base + std::string(config));
	rc->ReadConfig (cf.c_str());

	rc->Attribute ("maxgrad", &(p.mgr));
	rc->Attribute ("maxslew", &(p.msr));
	int ss ;
	rc->Attribute ("shots",   &ss);
	p.shots = ss;
	rc->Attribute ("dt",      &(p.dt));
	rc->Attribute ("res",     &(p.res));

	std::string rad (rc->GetText("/config/rad"));
	std::string fov (rc->GetText("/config/fov"));
	std::vector<std::string> rads = Split (rad, ",");
	std::vector<std::string> fovs = Split (fov, ",");

   	size_t nr = rads.size();
	assert (nr == fovs.size());

	p.fov = Matrix<double> (nr,1);
	p.rad = Matrix<double> (nr,1);
	for (size_t i = 0; i < nr; i++) {
		p.fov[i] = atof(fovs[i].c_str()); 
		p.rad[i] = atof(rads[i].c_str()); 
	}

	printf ("Computing variable density spiral for ... \n");
	printf ("    [maxgrad: %.2f, maxslew: %.2f, res: %.2f, dt: %.2e, shots: %li]\n", p.mgr, p.msr, p.res, p.dt, p.shots);
	p.res = p.res*.99;
	ticks start = getticks();
	Solution s = VDSpiral (p);
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
	


}
