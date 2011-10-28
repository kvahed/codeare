#include "GPUSimulator.hpp"

using namespace RRStrategy;

void 
GPUSimulator::SimulateExc  (const Matrix<cplx>&   txm, const Matrix<cplx>&   rf, 
							const Matrix<double>&  gr, const Matrix<double>&  r, 
							const Matrix<double>& b0m, const double&         dt, 
							const bool&             v, const size_t&        pos, 
							const int&            tid,       Matrix<double>&  m) {}


void 
GPUSimulator::SimulateRecv (const Matrix<cplx>&   rxm, const Matrix<double>& gr, 
							const Matrix<double>&   r, const Matrix<double>& m0, 
							const Matrix<double>& b0m, const double&         dt, 
							const bool&             v, const size_t&        pos, 
							const int&            tid,       Matrix<cplx>&  res) {}



void
GPUSimulator::Simulate     (const Matrix<cplx>&   txm, const Matrix<cplx>&   rxm, 
							const Matrix<cplx>&    rf, const Matrix<double>&  gr, 
							const Matrix<double>&   r, const Matrix<double>&  m0, 
							const Matrix<double>& b0m, const double&          dt, 
							const bool&           exc, const bool&             v, 
							const size_t&          np, 
						          Matrix<cplx>&   res, Matrix<double>&         m) {
	
	ticks            tic  = getticks();
	
	size_t           nr   =   r.Dim(1);
	size_t           ntxc = txm.Dim(1);
	size_t           nrxc = rxm.Dim(1);
	size_t           nt   =  gr.Dim(1);
	
	printf ("  Simulaing: %s on %04i isochromats ... ", (exc) ? "         excitation" : " signal acquisition", (int)nr); fflush(stdout);
	
	// Anything to do? ----------------
	
	if (gr.Size() < 1 || r.Size() < 1) {
		std::cout << "  Bailing out: %i Gradient step for %i isochromats? I don't think so!\n" << std::endl;
		return;
	}
	// --------------------------------
	
	Matrix<cplx> mres;
	
	if (!exc)
		mres = Matrix<cplx> (res.Dim(0), res.Dim(1), np);
	
	printf (" done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
	
}

