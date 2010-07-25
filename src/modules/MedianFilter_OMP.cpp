/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but 
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301  USA
 */

#include "MedianFilter_OMP.h"

#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_thread_num() { return 0;}
inline int omp_get_num_threads() { return 1;}
#endif

int compare_hints( const void* a, const void* b ) {
	int* arg1 = (int*) a;
	int* arg2 = (int*) b;
	if( *arg1 < *arg2 ) return -1;
	else if( *arg1 == *arg2 ) return 0;
	else return 1;
}

//process the median filter
RRSModule::error_code MedianFilter_OMP::ProcessData () {
	
	//default values
	const int m_ww = 25;  //window width
	const int m_wh = 25;  //window height
	
	int ex = (m_ww / 2), ey = (m_wh / 2);
	int array[m_ww*m_wh];
	
	cout << "MedianFilter::ProcessData() running on " << m_pixel.width() << "x" << 	m_pixel.height() << " image " << endl; 
	
	int x,y,fx,fy;
	
	//begin of parallel section: loop counters need private copies!
	
#pragma omp parallel default(shared)  private(y,fx,fy) 
	{
		
		int tid      = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int chunk    = m_pixel.width()/nthreads; //chunk-size for loop splitting
		
		if (tid==0) printf("\nNumber of threads = %d\n",nthreads);
		
#pragma omp for schedule(dynamic,chunk) 
		
		for (x=ex; x<m_pixel.width()-ex; ++x)
			for (y=ey; y<m_pixel.height()-ey; ++y) {
				for (fx=0; fx<m_ww; ++fx)
					for (fy=0; fy<m_wh; ++fy)
						array[fx+fy*m_ww] = m_pixel.at(x+fx-ex,y+fy-ey);
				qsort(array, m_ww*m_wh, sizeof(int), compare_hints );
				m_pixel.at(x,y) = array[m_ww*m_wh/2] ;
			}
		
	}
	
	return OK;

}

// the class factories
extern "C" ReconStrategy* create() {
    return new MedianFilter_OMP;
}

extern "C" void destroy (ReconStrategy* p) {
    delete p;
}
