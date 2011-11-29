void Rotate (const float* n, float* r)  {
	
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




void SimulateAcq (__global const float* rxm, __global const float*  gr, __global const float*   r, 
				  __global const float*   m, __global const float*  b0,          const float  gdt, 
				           const uint   pos, __global const uint*    n, __global       float* sig) {

	uint nt = n[0]; /* time points */
	uint nc = n[1]; /* channels    */

	float TWOPI = 6.283185307;

	float nv [3];    /* Rotation normal         */
	float rm [3*3];  /* Rotation matrix         */
	float lm [3];    /* Local magnetisation     */
	float tmp[3];    /* Temporary magnetisation */
	float lr [3];    /* Local position vector   */
	float ls [8][2]; /* Local sensitivities     */

    for (uint c = 0; c < nc; c++)  {
        ls[c][0] = rxm[2*pos*nc + c];
        ls[c][1] = rxm[2*pos*nc + c + 1];
	}
    
    for (size_t i = 0; i <  3; i++) {
        lr[i] = r[i + pos*3];
        lm[i] = m[i + pos*3]; 
    }

	nv[0] = 0.0;
	nv[1] = 0.0;

    for (size_t t = 0; t < nt; t++) {

        /* Acquisition: only gradients! */
        nv[2] = gdt * (gr[3*t+0]*lr[0] + gr[3*t+1]*lr[1] + gr[3*t+2]*lr[2] + b0[pos]*TWOPI);
        
		/* Calculate rotation */
        Rotate (nv, rm);
        
		/* Update m */
        tmp[0] = rm[0]*lm[0] + rm[3]*lm[1] + rm[6]*lm[2];
        tmp[1] = rm[1]*lm[0] + rm[4]*lm[1] + rm[7]*lm[2];
        tmp[2] = rm[2]*lm[0] + rm[5]*lm[1] + rm[8]*lm[2];
        
        lm [0] = tmp[0];
        lm [1] = tmp[1];
        lm [2] = tmp[2];
        
        /* Weighted signal */
        for (size_t c = 0; c < nc; c++) {
            sig[t + 2*nc*c    ] += ls[c][0]*lm[0];
			sig[t + 2*nc*c + 1] += ls[c][1]*lm[1];
		}
        
    }

}




void SimulateExc  (__global const float* txm, __global const float*  gr, __global const float* rf, 
				   __global const float*   r, __global const float*  b0, __global const float gdt, 
				            const uint   pos,          const uint*    n, __global       float*   m) {

	uint nt = n[0]; /* time points */
	uint nc = n[1]; /* channels    */

	float TWOPI = 6.283185307;

    float nv [3];    /* Rotation axis        */
	float rm [3*3];  /* Rotation matrix      */
    float lm [3];    /* Local magnetisation  */
    float tmp[3];    /* Temp magnetisation   */
    float lr [3];    /* Local spatial vector */
    float ls [8][2]; /* Local sensitivities  */
	float rfs[2];    /* Local composite RF   */

    for (uint i = 0; i < 3; i++) lr[i] = r[i+pos*3];

    for (uint i = 0; i < nc; i++) {
		ls[i][0] = txm[2*pos*nc + i];
        ls[i][1] = txm[2*pos*nc + i + 1];
	}

	/* Equilibrium */
	lm[0] = 0.0;
	lm[1] = 0.0;
    lm[2] = 1.0;
    
    for (size_t t = 0; t < nt; t++) {
        
        size_t rt = nt-1-t;
        
        rfs[0] = 0.0;
		rfs[1] = 0.0;
        
        for (size_t i = 0; i < nc; i++) {
			rfs[0] += (rf[rt + i*nt*2  ] * ls[i][0] - rf[rt + i*nt*2+1] * ls[i][1]);
			rfs[1] += (rf[rt + i*nt*2+1] * ls[i][0] + rf[rt + i*nt*2  ] * ls[i][1]);
		}

        nv[0] = gdt * -rfs[0];
        nv[1] = gdt *  rfs[1];
        nv[2] = gdt * (- gr[rt*3+0] * lr[0] - gr[rt*3+1] * lr[1] - gr[rt*3+2] * lr[2] + b0[pos]*TWOPI);
        
        Rotate (nv, rm);
        
        tmp[0] = rm[0]*lm[0] + rm[3]*lm[1] + rm[6]*lm[2];
        tmp[1] = rm[1]*lm[0] + rm[4]*lm[1] + rm[7]*lm[2];
        tmp[2] = rm[2]*lm[0] + rm[5]*lm[1] + rm[8]*lm[2];
        
        lm[0]  = tmp[0];
        lm[1]  = tmp[1];
        lm[2]  = tmp[2];
        
    }
    
    m[pos*3  ] = lm[0];
    m[pos*3+1] = lm[1];
    m[pos*3+2] = lm[2]; 

}



void CollectSignals (__global T *g_idata, __global T *g_odata, unsigned int n, __local volatile T* sdata) {


    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = get_local_id(0);
    unsigned int i   = get_group_id(0)*(get_local_size(0)*2) + get_local_id(0);

    sdata[tid] = (i < n) ? g_idata[i] : 0;
    if (i + get_local_size(0) < n) 
        sdata[tid] += g_idata[i+get_local_size(0)];  

    barrier(CLK_LOCAL_MEM_FENCE);

    // do reduction in shared mem
    #pragma unroll 1
    for(unsigned int s=get_local_size(0)/2; s>32; s>>=1) 
    {
        if (tid < s) 
        {
            sdata[tid] += sdata[tid + s];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (tid < 32)
    {
        if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; }
        if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; }
        if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; }
        if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; }
        if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; }
        if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; }
    }

    // write result for this block to global mem 
    if (tid == 0) g_odata[get_group_id(0)] = sdata[0];

}



__kernel void Simulate     (__global const float*  b1, __global const float*  gr, __global const float*  sr, 
							__global const float*  b0, __global const float*  tm, __global const float*  sm, 
							__global const float* jac, __global const float* gdt, __global const  uint*   n,
							__global       float*  rf, __global       float*   m) 
{

	uint igl = get_global_id(0);
	uint igr = get_group_id (0);
	uint ilc = get_local_id (0);

	uint nt = n[0]; /* # timepoints  */
	uint nc = n[1]; /* # channels    */
	uint nr = n[2]; /* # Sites   */

	__global float sig[2 * nt * nc];

	SimulateAcq (b1, gr, r, tm, b0, gdt[0], igl, n, sig);
	barrier(CLK_LOCAL_MEM_FENCE);

	CollectSignals();
	barrier(CLK_LOCAL_MEM_FENCE);

	SimulateExc (b1, gr, rf, r, b0, gdt[0], igl, n, m);
	barrier(CLK_LOCAL_MEM_FENCE);

}

