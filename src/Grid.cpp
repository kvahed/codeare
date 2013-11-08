#include "Grid.hpp"
#include "BLACS.hpp"
#include <sstream>
#include <math.h>
#include <iostream>

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

static inline void 
grid_init (Grid& gd) {

	/* General sizes */
	Cblacs_pinfo(&gd.rk, &gd.np);

	/* Grid setup */
	grid_setup (gd.np, gd.nr, gd.nc);
	Cblacs_get (-1, 0, &gd.ct);
	Cblacs_gridinit (&gd.ct, &gd.order, gd.nr, gd.nc);
    Cblacs_gridinfo (gd.ct, &gd.nr, &gd.nc, &gd.mr, &gd.mc);

}


static inline void 
grid_exit (Grid& gd) {

	/* Clean up */
	Cblacs_gridexit (gd.ct);
	Cblacs_exit (0);

}
#endif 


Grid* Grid::m_inst = 0;


Grid::Grid () :
    np(0), rk(0), ct(0), nr(0), nc(0), mr(0), mc(0), order('R') {
    
#ifdef HAVE_MPI
    grid_init (*this);
    
#ifdef BLACS_DEBUG2
    std::cout << *this;
#endif
    
#endif
    
}


Grid::~Grid () {
    
#ifdef HAVE_MPI
    grid_exit(*this);
#endif
    
}

const char* Grid::c_str() const {
    
    std::stringstream ss;
    ss << "- np(" << np << ") rk(" << rk << ") ct(" << ct << ") nr(" << nr
       << ") nc(" << nc << ") mr(" << mr << ") mc(" << mc << ") or(" << order
       << ")\n";
    return ss.str().c_str();
    
}

Grid& Grid::Instance() {
    
    
    if (m_inst == 0)
        m_inst = new Grid ();
    
	return *m_inst;

}

