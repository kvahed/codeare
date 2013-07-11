#include "../modules/VDSpiral.hpp"
#include "Statistics.hpp"
#include "Tokenizer.hpp"
#include "SimpleTimer.hpp"

#include <sstream>

template <class T> bool
vdspiraltest (Connector<T>* rc) {

	SpiralParams p;

	rc->Attribute ("maxgrad", &(p.mgr));
	rc->Attribute ("maxslew", &(p.msr));
	int ss ;
	rc->Attribute ("shots",   &ss);
	p.shots = ss;
	rc->Attribute ("dt",      &(p.dt));
	rc->Attribute ("res",     &(p.res));

	rc->Attribute ("gunits",  &(p.gunits));
	rc->Attribute ("lunits",  &(p.lunits));

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

	printf ("\nComputing variable density spiral for ... \n");
	printf ("    [maxgrad: %.2f, maxslew: %.2f, res: %.2f, dt: %.2e, shots: %li]\n\n",
            p.mgr, p.msr, p.res, p.dt, p.shots);

	ticks start = getticks();
    SimpleTimer st ("VD spiral design");
	Solution s = VDSpiral (p);
    st.Stop();

    std::stringstream ofname;
    ofname << base;
    ofname << "/";
    ofname << rc->GetElement("/config/data-out")->Attribute("fname");
    
    IOContext f = fopen (ofname.str(), WRITE);
    s.dump (f);
    fclose (f);
	
	return true;

}
