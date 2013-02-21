#include "Workspace.hpp"
#include "BLACS.hpp"

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




Workspace* Workspace::m_inst = 0; 


Workspace::Workspace() {

    // Defaults
    m_gd.rk = 0; 
    m_gd.np = 0; 
    m_gd.mr = 0; 
    m_gd.mc = 0; 
    m_gd.nc = 0; 
    m_gd.nr = 0; 
    m_gd.ct = 0;
    
#ifdef HAVE_MPI
    grid_init (m_gd);
    
#ifdef BLACS_DEBUG2
    std::cout << m_gd;
#endif
#endif
    
}


Workspace::~Workspace () { 
    
	Finalise();
	m_inst = 0;
    
#ifdef HAVE_MPI
    grid_exit(m_gd);
#endif
	
}


Workspace* 
Workspace::Instance ()  {

    if (m_inst == 0) 
        m_inst = new Workspace ();

	bool t = true;
	bool f = false;
	
	m_inst->SetAttribute("FFTWFThreadsInitialised", &f);
	m_inst->SetAttribute("FFTWThreadsInitialised", &f);
	
    return m_inst;
		
}
	
error_code
Workspace::Finalise () {
	
	while (!m_ref.empty()) {

  		reflist::iterator nit = m_ref.begin();
  		store::iterator dit = m_store.find (nit->second[0]);

		if      (nit->second[1].compare(typeid(cxfl).name())   == 0)
			delete boost::any_cast<Ptr<Matrix<cxfl  > > >(dit->second);
		else if (nit->second[1].compare(typeid(cxdb).name())   == 0)
			delete boost::any_cast<Ptr<Matrix<cxdb  > > >(dit->second);
		else if (nit->second[1].compare(typeid(float).name())  == 0)
			delete boost::any_cast<Ptr<Matrix<float > > >(dit->second);
		else if (nit->second[1].compare(typeid(double).name()) == 0)
			delete boost::any_cast<Ptr<Matrix<double> > >(dit->second);
		else if (nit->second[1].compare(typeid(short).name())   == 0)
			delete boost::any_cast<Ptr<Matrix<short > > >(dit->second);
		else if (nit->second[1].compare(typeid(long).name())   == 0)
			delete boost::any_cast<Ptr<Matrix<long  > > >(dit->second);

  		m_store.erase(dit);
  		m_ref.erase(nit);

  	}

	return OK;
	
}

