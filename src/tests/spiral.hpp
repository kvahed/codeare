#include "../modules/VDSpiral.hpp"

template <class T> bool
vdspiraltest (Connector<T>* rc) {

	SpiralParams p;

	Matrix<double> fov (4,1);
	Matrix<double> rad (4,1);
	
	fov[0] = 35.0; rad[0] = 0.0;
	fov[1] = 35.0; rad[1] = 0.1;
	fov[2] = 10.0; rad[2] = 0.13;
	fov[3] = 10.0; rad[3] = 1.0;

	p.fov = &fov;
	p.rad = &rad;
			
	p.mgr =  4.0;
	p.msr = 20.0;
	p.dt  = 1.0e-2;
	p.shots = 1;
	p.res = 6.0;
	
	Spiral s = VDSpiral (p);


}
