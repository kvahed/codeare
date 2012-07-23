#include <math.h>
#include "PMatrix.hpp"
#include "PIO.hpp"

#ifdef HAVE_MPI
static void 
grid_setup (const int& np, int& nc, int& nr) {
	
	int sqrtnp, i;
	
	sqrtnp = int(sqrt(double(np)) + 1);
	
	for (i = 1; i < sqrtnp; i++)
		if (!(np % i))
			nc = i;
	
	nr = np/nc;
	
}


static void 
sl_init (grid_dims& gd) {

	/* General sizes */
	Cblacs_pinfo(&gd.rk, &gd.np);

	/* Grid setup */
	grid_setup (gd.np, gd.nr, gd.nc);
	Cblacs_get (-1, 0, &gd.ct);
	Cblacs_gridinit (&gd.ct, &gd.order, gd.nr, gd.nc);

}


static void 
sl_finalise (grid_dims& gd) {

	/* Clean up */
	Cblacs_gridexit (gd.ct);
	Cblacs_exit (0);

}
#endif 


template <class T> bool
mpitest (Connector<T>* rc) {

#ifdef HAVE_MPI
	grid_dims gd;
	gd.order = 'C';

	/* Initialise ScaLAPACK */
	sl_init (gd);

	/* MPI aware matrix */
	PMatrix<double> pm (32,40);
	pm.Val() = 20.0;
	print (pm);

	/* Finalise ScaLAPACK */
	sl_finalise (gd);
#else 
	printf ("MPI not supported by this build\n");
#endif

	return 0;

}
