#include "../modules/VDSpiral.hpp"

template <class T> bool
vdspiraltest (Connector<T>* rc) {

	SpiralParams p;

	/*
	p.fov = Matrix<double> (4,1);
	p.rad = Matrix<double> (4,1);
	p.fov[0] = 35.0; p.rad[0] = 0.0;
	p.fov[1] = 35.0; p.rad[1] = 0.1;
	p.fov[2] = 10.0; p.rad[2] = 0.13;
	p.fov[3] = 10.0; p.rad[3] = 1.0;
	*/

	p.fov = Matrix<double> (2,1);
	p.rad = Matrix<double> (2,1);
	p.fov[0] = 19.2; p.rad[0] = 0.0;
	p.fov[1] = 19.2; p.rad[1] = 1.0;

	p.mgr =   4.0;
	p.msr =  20.0;
	p.dt  =   1.0e-2;
	p.shots = 3;
	p.res =   6.0;
	
	printf ("Computing variable density spiral ...\n");
	ticks start = getticks();
	Spiral s = VDSpiral (p);
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
