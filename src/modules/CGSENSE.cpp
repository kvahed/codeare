#include "CGSENSE.hpp"
#include "nfftstub.h"

#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_thread_num() { return 0;}
inline int omp_get_num_threads() { return 1;}
#endif

#include <vector>
#include <cycle.h>

#include <stdint.h>
#include <sys/types.h>
#include <sys/sysctl.h>

#define M 714025
#define IA 1366
#define IC 150889
#define LM 2147483647
#define LAM (1.0/LM)
#define LA 16807
#define LR 2836
#define LQ 127773

long idum;

std::string sides[3] = {"Nx", "Ny", "Nz"};

using namespace RRStrategy;

CGSENSE::CGSENSE () {

}

RRSModule::error_code 
CGSENSE::Init() {
	
	RRSModule::error_code error = OK; 

	m_N   = new int[3];
	m_n   = new int[3];

	for (int i = 0; i < 3; i++) {
		m_N[i] = 0; 
		m_n[i] = 0;
	}

	m_testcase = 0;
	m_verbose  = 0;
	m_noise    = 0;

	// Verbosity ----------------------------

	Attribute ("testcase",  &m_testcase);
	// --------------------------------------

	// Verbosity ----------------------------

	Attribute ("verbose",   &m_verbose);
	// --------------------------------------

	// Noise --------------------------------

	Attribute ("noise",   &m_noise);
	// --------------------------------------

	// Dimensions ---------------------------

	Attribute("dim",       &m_dim);

	for (int i = 0; i < m_dim; i++)
		Attribute (sides[i].c_str(),       &m_N[i]);

	m_M   = m_helper.Size();
	// --------------------------------------

	// iNFFT convergence and break criteria -

	Attribute ("maxit",   &m_maxit);
	Attribute ("epsilon", &m_epsilon);
	// --------------------------------------

	// CG convergence and break criteria

	Attribute ("cgeps",   &m_cgeps);
	Attribute ("cgmaxit", &m_cgmaxit);
	// --------------------------------------

	// Oversampling -------------------------

	int      m           = 0;
	double   alpha       = 0.0;

	Attribute ("m",       &m);
	Attribute ("alpha",   &alpha);

	for (int i = 0; i < m_dim; i++)
		m_n[i] = ceil (m_N[i]*alpha);
	// --------------------------------------

	// Initialise FT plans ------------------
	
	for (int i = 0; i < NTHREADS; i++)
		nfft::init (m_dim, m_N, m_M, m_n, m, &m_fplan[i], &m_iplan[i], m_epsilon);
	// --------------------------------------

	// Allocate RAM for fix size memories (weights and k-space)

	m_ftk      = (double*) malloc (    m_dim  * m_M    * sizeof(double)); 
	m_ftw      = (double*) malloc (             m_M    * sizeof(double)); 
	// --------------------------------------

	return error;

}


/**
 * @brief               Compute left hand side (i.e. Multiply E with spatial (image) data)
 *                      Forward NFFT in and elementwise multiply with spatial sensitivity of every channel 
 * 
 * @param  in           Original discretised sample O (Nx x Ny x Nz)
 * @param  sm           Sensitivity maps            O (Nx x Ny x Nz x Nc)
 * @param  np           Non-Cartesian strategy for non uniform ft
 * @param  out          Result                      O (Nk x Nc)
 */
RRSModule::error_code 
E  (Matrix<raw>* in, Matrix<raw>* sm, nfft_plan* np, Matrix<raw>* out, int dim) {

	// Clear output container
	out->Zero();
	
	// Some dimensions
	int        ncoils   = sm->Dim(dim);
	int        nsamples = out->Size() / ncoils;
	int        imgsize  = in->Size();

	// Create container for FT input
	double*    ftout    = (double*) malloc (2 * ncoils * nsamples * sizeof(double));

	// Loop over coils, Elementwise multiplication of maps with in (s.*in), ft and store in out
	
#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		omp_set_num_threads(NTHREADS);
		
#pragma omp for
		for (int j = 0; j < ncoils; j++) {
			
			double*    ftin     = (double*) malloc (2 * imgsize  * sizeof(double));
			
			int     ipos     = j * imgsize;
			int     spos     = j * nsamples;
			
			// Copy data to FT
			for (int i = 0; i < imgsize; i++) {
				
				raw tmp     = sm->at(ipos + i) * in->at(i);
				
				ftin[2*i  ] = tmp.real(); 
				ftin[2*i+1] = tmp.imag(); 
				
			}
			
			// Forward ft
			nfft::ft (&np[tid], ftin, &ftout[2*spos]);
			
			free (ftin);
			
			// Copy FTed data back
			for (int i = 0; i < nsamples; i++) 
				out->at(spos + i) = raw(ftout[2*i+2*spos],ftout[2*i+1+2*spos]);
			
		}
	}
	
	// Free RAM
	free (ftout);
	
	// Return success
	return OK;
	
}

struct CpuInfo {
	char vendor_id[50];
	int family;
	char model[50];
	float freq;
	char cache[20];
};

double clockrate () {
#if defined(HAVE_MACH_ABSOLUTE_TIME)

    uint64_t freq = 0;
    size_t   size = sizeof(freq);
	
    if (sysctlbyname("hw.tbfrequency", &freq, &size, NULL, 0) < 0)
		perror("sysctl");
	
    return freq;

#else

	struct CpuInfo info = {"", 0, "", 0.0, ""};
	
	FILE *cpuInfo;
	
	std::string fname = "/proc/cpuinfo";
	
	if ( ( cpuInfo = fopen(fname.c_str(), "rb") == NULL ) )
		printf("Error! Cannot open %s", fname.c_str());
	
	else 
		while (!feof(cpuInfo)) {

			fread(&info, sizeof(struct CpuInfo), 1, cpuInfo);
			if(info.family !=0) {
				printf("%s\n%d\n%s\n%.2f\n%s\n", info.vendor_id, info.family, info.model, info.freq, info.cache);
			}
		}

	return (double)info.freq;

#endif
}

/**
 * @brief               Compute right hand side (i.e. Multiply E^H, Hermitian counterpart to E, with k-space data)
 *
 * @param  in           K-space samples along trajectory O (Nk x Nc)
 * @param  sm           Sensitivity maps                 O (Nx x Ny x Nz x Nc)
 * @param  np           NuFFT plan
 * @param  spc          Solver plan
 * @param  epsilon      Convergence criterium for ift (default 3e-7)
 * @param  maxit        Maximum number of solver iterations (default 3)
 * @param  out          Returned product                 O (Nx x Ny x Nz)
 */
RRSModule::error_code
EH (Matrix<raw>* in, Matrix<raw>* sm, nfft_plan* np, solver_plan_complex* spc, double epsilon, int maxit, Matrix<raw>* out, int dim) {

	// Clear outgoing container
	out->Zero();

	// Some dimensions
	int        ncoils   = sm->Dim(dim);
	int        nsamples = in->Size() / ncoils;
	int        imgsize  = out->Size();

	double* ftout = (double*) malloc (2 * ncoils * imgsize  * sizeof(double));
	
	// OMP Loop over coils, Inverse FT every signal in *in, 
	// Sum elementwise mutiplied images with according sensitivity maps 
#pragma omp parallel default (shared) 
	{
		
		int tid      = omp_get_thread_num();
		omp_set_num_threads(NTHREADS);
		
#pragma omp for
		for (int j = 0; j < ncoils; j++) {
			
			// Containers for FT I/O
			double* ftin  = (double*) malloc (2 * nsamples * sizeof(double));
			
			int    spos   = j * nsamples;
			int    ipos   = j * imgsize;
			
			// Copy to iFT
			for (int i = 0; i < nsamples; i++) {
				ftin[2*i  ] = (in->at(spos + i)).real();
				ftin[2*i+1] = (in->at(spos + i)).imag();
			}
			
			// Inverse FT
			nfft::ift (&np[tid], &spc[tid], ftin, &ftout[2*ipos], maxit, epsilon);
			
			free (ftin);
			
		}
	}

	for (int j = 0; j < ncoils; j++) {
		int    ipos   = j * imgsize;
		for (int i = 0; i < out->Size(); i++) {
			raw sens = sm->at(ipos + i);
			out->at(i) += raw(ftout[2*i+2*ipos], ftout[2*i+1+2*ipos]) * conj(sens);
		}
	}
	
	// Free RAM
	free (ftout);

	return OK;
	
}

//--------------------------------------------------------------------
//       Returns uniform deviate between 0.0 and 1.0.
//       Used to generate PN data
//---------------------------------------------------------------------
void setseed (double i=1349555.0) { idum = (long)i;}

double uniform () {
	
	double  r1;
	long    hi;
	hi      = idum/LQ;
	idum    = LA*(idum-hi*LQ) - LR*hi;
	
	if (idum < 0) 
		idum += LM;
	
	r1  = LAM*idum;

	return(r1);

}

raw WhiteNoise () {

	float fac, r, v1, v2;

	do {
		v1 = (2.0 * uniform()) - 1.0;
		v2 = (2.0 * uniform()) - 1.0;
		r = (v1*v1) + (v2*v2);
	} while (r >= 1.0);

	fac = sqrt(-2.0 * log(r) / r);

	return ( raw( v2*fac, v1*fac ) );
}

RRSModule::error_code
AddPseudoRandomNoise (Matrix<raw>* m, float max) {

	setseed();

	for (int i = 0; i < m->Size(); i++)
		m->at(i) = m->at(i) + max*WhiteNoise();
	
	return OK;

}


RRSModule::error_code
CGSENSE::Process () {

	RRSModule::error_code error = OK;


	Matrix<raw> a, p, q, r, r_new;

	if (m_noise > 0.0)
		AddPseudoRandomNoise (&m_raw, (float)m_noise);

	for (int i = 0; i < INVALID_DIM; i++)
		a.Dim(i) = 1;
	for (int i = 0; i < m_dim      ; i++)
		a.Dim(i) = m_N[i];

	a.Reset();

	memcpy (m_ftw, &m_helper[0], m_helper.Size()*sizeof(double));
	memcpy (m_ftk, &m_kspace[0], m_kspace.Size()*sizeof(double));

	for (int i = 0; i < NTHREADS; i++) {
		nfft::kspace  (&m_fplan[i],              m_ftk);
		nfft::weights (&m_fplan[i], &m_iplan[i], m_ftw);
	}

	m_sens = m_rhelper;

	Matrix<raw> store;
	for (int i = 0; i < INVALID_DIM; i++)
		store.Dim(i) = 1;

	if (m_verbose == 1) 
		for (int i = 0; i < m_dim; i++)
			store.Dim(i) = m_N[i];

	store.Dim(m_dim) = m_cgmaxit;
	store.Reset();

	m_rhelper = m_raw;

	Matrix<raw> sigtmp;
	sigtmp.Dim (COL) = m_M;
	sigtmp.Dim (LIN) = 8;
	sigtmp.Reset();

	Matrix<raw> imgtmp;
	
	for (int i = 0; i < m_dim; i++)
		imgtmp.Dim (i) = m_N[i];

	imgtmp.Reset();

	if (m_testcase) {
		E  (&m_rhelper, &m_sens, &m_fplan[0], &sigtmp, m_dim);
		m_rhelper = sigtmp;
		m_raw     = sigtmp;
	}

	ticks cgstart = getticks();

	EH (&m_raw, &m_sens, &m_fplan[0], &m_iplan[0], m_epsilon, m_maxit, &a, m_dim);
	p = a;
	r = a;
	q = a;

	// Out going image
	// Resize m_raw for output
	for (int i = 0; i < INVALID_DIM; i++)
		m_raw.Dim(i) = 1;

	for (int i = 0; i < m_dim; i++)
		m_raw.Dim(i) = m_N[i];

	m_raw.Reset();

	std::vector<double> res;
	
	float       rn    = 0.0;
	float       an    = 0.0;
	float       rnewn = 0.0;
	raw         rtmp  = raw(0.0,0.0);
	int         iters = 0;

	// CG iterations (Pruessmann et al. (2001). MRM, 46(4), 638-51.)
	for (int i = 0; i < m_cgmaxit; i++, iters++) {

		rn = r.norm().real();
		an = a.norm().real();
		
		res.push_back(rn/an);
		
		printf ("%03i: CG residuum: %.9f\n", i, res.at(i));
		if (res.at(i) <= m_cgeps)
			break;
		
		E  (&p,      &m_sens, &m_fplan[0],                               &sigtmp, m_dim);
		EH (&sigtmp, &m_sens, &m_fplan[0], &m_iplan[0], m_epsilon, m_maxit, &q     , m_dim);
		
		rtmp      = (rn / (p.dotc(q)));
		imgtmp    = p * rtmp;
		m_raw     = m_raw + imgtmp;
		sigtmp    = sigtmp * rtmp;
		m_rhelper = m_rhelper + sigtmp;
		imgtmp    = q * rtmp;
		r_new     = r - imgtmp;
		rnewn     = r_new.norm().real();
		rtmp      = rnewn/rn;
		imgtmp    = p * rtmp;
		p         = r_new + imgtmp;
		r         = r_new;

		if (m_verbose == 1)
			memcpy (&store[i*m_raw.Size()], &m_raw[0], m_raw.Size() * sizeof(double));

	}

	printf ("Processing CG-SENSE took: %.4f seconds.\n", elapsed(getticks(), cgstart) / clockrate());

	if (m_verbose == 1) {
		m_raw.Dim(2) = iters;
		m_raw.Reset();
		memcpy (&m_raw[0], &store[0], m_raw.Size() * sizeof(double));
	}

	return error;

}


// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new CGSENSE;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}

