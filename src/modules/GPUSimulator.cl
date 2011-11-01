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


void SimulateAcq (const float* rxm, const float*  gr, const float* r, 
				  const float*   m, const float* b0m, const float* dt, 
				  const int*   pos,       float* sig) {

	uint nd = 3;

	uint nt = n[0]; /* time points */
	uint nc = n[1]; /* channels    */

	float nv [nd];
	float rm [nd*nd];
	float ml [nd];
	float tmp[nd];
	float lr [nd];
	float ls [nc][2];

    for (uint c = 0; c < nc; c++)  {
        ls[c][0] = rxm[2*pos*nc + c];
        ls[c][1] = rxm[2*pos*nc + c + 1];
	}
    
    for (size_t i = 0; i <  nd; i++) {
        lr[i] = r[i + pos*nd];
        lm[i] = m[i + pos*nd]; 
    }

	nv[0] = 0.0;
	nv[1] = 0.0;

    for (size_t t = 0; t < nt; t++) {

        /* Acquisition: only gradients! */
        nv[2] = gdt * (gr[nd*t+0]*lr[0] + gr[nd*t+1]*lr[1] + gr[nd*t+2]*lr[2] + b0[pos]*TWOPI);
        
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
			sig[t + 2*nc*c + 1] += ls[c][1]*lm[1]
		}
        
    }

}

void SimulateExc  (const float* txm, const float*  gr, const float* rf, 
				   const float*   r, const float* b0m, const float* dt, 
				   const int*   pos,       float*   m) {

	uint nt = n[0]; /* time points */
	uint nc = n[1]; /* channels    */
	uint nr = n[3]; /* Exc sites   */
	uint nd = 3;

    float nv [nd];    /* Rotation axis        */
	float rm [nd*nd]; /* Rotation matrix      */
    float ml [nd];    /* Local magnetisation  */
    float tmp[3];     /* Temp magnetisation   */
    float lr [3];     /* Local spatial vector */
    float ls [nc][2]; /* Local sensitivities  */
	float rfs[2];     /* Local composite RF   */

    for (uint i = 0; i < nd; i++) lr[i] = r  (i,pos);

    for (uint i = 0; i < nc; i++) {
		ls[c][0] = txm[2*pos*nc + c];
        ls[c][1] = txm[2*pos*nc + c + 1];
	}

	/* Equilibrium */
	ml[0] = 0.0;
	ml[1] = 0.0;
    ml[2] = 1.0;
    
    for (size_t t = 0; t < nt; t++) {
        
        size_t rt = nt-1-t;
        
        rfs[0] = 0.0;
		rfs[1] = 0.0;
        
        for (size_t i = 0; i < nc; i++) {
			rfs[0] += (rf[rt + i*nt*2  ] * ls[c][0] - rf[rt + i*nt*2+1] * ls[c][1]);
			rfs[1] += (rf[rt + i*nt*2+1] * ls[c][0] + rf[rt + i*nt*2  ] * ls[c][1]);
		}

        n[0] = gdt * -rfs[0];
        n[1] = gdt *  rfs[1];
        n[2] = gdt * (- gr(rt*nd+0) * lr[0] - gr(rt*nd+1) * lr[1] - gr(rt*nd+2) * lr[2] + b0[pos]*TWOPI);
        
        Rotate (nv, rm);
        
        tmp[0] = rm[0]*ml[0] + rm[3]*ml[1] + rm[6]*ml[2];
        tmp[1] = rm[1]*ml[0] + rm[4]*ml[1] + rm[7]*ml[2];
        tmp[2] = rm[2]*ml[0] + rm[5]*ml[1] + rm[8]*ml[2];
        
        ml[X]  = tmp[0];
        ml[Y]  = tmp[1];
        ml[Z]  = tmp[2];
        
    }
    
    m[pos*nd  ] = ml[0];
    m[pos*nd+1] = ml[1];
    m[pos*nd+2] = ml[2]; 

}



void ReduceSignals (const float* data, volatile float* rdata) {
}


__kernel void Simulate     (__global const float* rxm, __global const    float* txm, 
							__global const float*  gr, __global const    float*  tr, 
							__global const float*  sr, __global const    float* tb0, 
							__global const float* sb0, __global const    float*  tm,
							__global const float*  sm, __global const    float* jac, 
							__global const float* gdt, __global const      int*   n,
							__global const int*    nc, __global const      int*  nr,
							__global const int*    
							__global       float*  rf, __global          float*   m) 
{

	uint igl = get_global_id(0);
	uint igr = get_group_id (0);
	uint ilc = get_local_id (0);

	__local float[nt*nc] sig;

	SimulateAcq (rxm, gr, tr, tm, tb0, gdt, igl, n, sig);
	
	barrier(CLK_LOCAL_MEM_FENCE);

	ReduceSignals();

	barrier(CLK_LOCAL_MEM_FENCE); 

	SimulateExc (txm, gr, rf, sr, sb0, gdt, igl, n, m);

}

