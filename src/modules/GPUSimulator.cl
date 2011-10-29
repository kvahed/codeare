#define STRINGIFY(A) #A

std::string kernel_source = STRINGIFY
(
									  
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
