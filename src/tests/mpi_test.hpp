#include <math.h>
#include "PMatrix.hpp"
#include "PIO.hpp"

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


int main (int argc, char** argv) {

	grid_dims gd;
	gd.order = 'C';

	/* Initialise ScaLAPACK */
	sl_init (gd);

	PMatrix<double> pm (32,40);
	pm.Val() = 20.0;
	print (pm);

	/* Finalise ScaLAPACK */
	sl_finalise (gd);

	return 0;

}
