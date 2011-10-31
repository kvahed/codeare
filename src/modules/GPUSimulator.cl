void Rotate (const float* n, float* r) 
{
	
	float phi = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
	
	if (!phi) {
		
		r [0] = 1.0;
		r [1] = 0.0;
		r [2] = 0.0;
		r [3] = 0.0;
		r [4] = 1.0;
		r [5] = 0.0;
		r [6] = 0.0;
		r [7] = 0.0;
		r [8] = 1.0;
		
	} else {
		
		float hp    =  phi/2;		
		float cp    =  cos(hp);
		float sp    =  sin(hp)/phi;
		float ar    =  cp;
		float ai    = -n[2]*sp;
		float br    =  n[1]*sp;
		float bi    = -n[0]*sp;
		
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


void SimulateRecv (const float* rxm, const float* gr, const float* r, const float*  m0, 
			  const float* b0m, const float* dt, const int* pos,       float* sig) 
{
}

void SimulateExc  (const float* txm, const float* rf, const float* gr, const float*  r, 
			  const float* b0m, const float* dt, const int*  pos,       float*  m) 
{
}

__kernel void Simulate     (__global const float* rxm, __global const float* txm, 
							__global const float*  gr, __global const float*  tr, 
							__global const float*  sr, __global const float* tb0, 
							__global const float* sb0, __global const float*  tm,
							__global const float*  sm, __global const float* jac, 
							__global const float* gdt, __global const   int*  nt,
							__global const int*    nc, 
							__global       float*  rf, __global       float*   m) 
{

    m_cmdq.finish();

}

