#define STRINGIFY(A) #A

std::string kernel_source = STRINGIFY
	(
	 
	 
	 void 
	 Rotate (const float* n, float* r) {
		 
		 float phi = sqrt(n[X]*n[X] + n[Y]*n[Y] + n[Z]*n[Z]);
		 
		 // Identity
		 if (!phi)

			 r = Matrix<float>::Id(3);
		 
		 else {
			 
			 // Cayley-Klein parameters
			 float hp    =  phi/2;		
			 float cp    =  cos(hp);
			 float sp    =  sin(hp)/phi; // n is unit length.
			 float ar    =  cp;
			 float ai    = -n[Z]*sp;
			 float br    =  n[Y]*sp;
			 float bi    = -n[X]*sp;
			 
			 float arar  =   ar*ar;
			 float aiai  =   ai*ai;
			 float arai2 = 2*ar*ai;
			 float brbr  =   br*br;
			 float bibi  =   bi*bi;
			 float brbi2 = 2*br*bi;
			 float arbi2 = 2*ar*bi;
			 float aibr2 = 2*ai*br;
			 float arbr2 = 2*ar*br;
			 float aibi2 = 2*ai*bi;
			 
			 r[0] =  arar  - aiai - brbr + bibi;
			 r[1] = -arai2 - brbi2;
			 r[2] = -arbr2 + aibi2;
			 r[3] =  arai2 - brbi2; 
			 r[4] =  arar  - aiai + brbr - bibi;
			 r[5] = -aibr2 - arbi2;
			 r[6] =  arbr2 + aibi2;
			 r[7] =  arbi2 - aibr2;
			 r[8] =  arar  + aiai - brbr - bibi;
			 
		 }
		 
	 }
	 
	 
	 void 
	 SimulateRecv (const float* rxm, const float* gr, const float* r, const float*  m0, 
				   const float* b0m, const float  dt, const int  pos,       float* sig) {}
	 
	 void 
	 SimulateExc  (const float* txm, const float* rf, const float* gr, const float*  r, 
				   const float* b0m, const float  dt, const int   pos,       float*  m) {}
	 
	 __kernel void 
	 Simulate     (__global const float* rxm, __global const float* txm, 
				   __global const float*  rf, __global const float*  gr,
				   __global const float*   r, __global const float*  m0,
				   __global const float* b0m, __global const float   dt,
				   __global const bool   exc, __global const float  res) {}
	 
	 );
