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

template <class T> T** CreateImage(int n,int m){
	T** pp = new T* [n];
	pp[0] = new T [n*m];
	for (int i = 1; i < n; ++i) pp[i] = pp[i-1] + m;
	return pp;
}

int t_compare_ints( const void* a, const void* b ) {
	int* arg1 = (int*) a;
	int* arg2 = (int*) b;
	if( *arg1 < *arg2 ) return -1;
	else if( *arg1 == *arg2 ) return 0;
	else return 1;
}

template <class T> void DeleteImage(T** pp){
	delete [] pp[0];
	delete [] pp;
}


RRSModule::error_code MedianFilter_OMP::Process () {
	
	const int ww = 25;
	const int wh = 25;

	int iw = m_pixel.width(); //image width 
	int ih = m_pixel.height(); //image height

	int ex = (ww / 2), ey = (wh / 2);
	int array[ww*wh];
	
	int x,y,fx,fy;
	int** input_image  = CreateImage<int>(iw,ih);


#pragma omp parallel default(shared)  private(y,fx,fy)
	{
		int tid      = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int chunk    = iw/nthreads; //chunk-size for loop splitting

		if (tid==0) 
			cout << "MedianFilter_OMP::Process() running on " << iw << "x" << ih << " image with " << nthreads << " threads" << endl; 

#pragma omp for schedule(dynamic,chunk)
		for (x=0; x<iw; ++x)
			for (y=0; y<ih; ++y)
					input_image[x][y]  = m_pixel.at(x,y);

		int** array = CreateImage<int>(ww,wh); //local to each thread !

#pragma omp for  schedule(dynamic,chunk)
		for (x=ex; x<iw-ex; ++x)
			for (y=ey; y<ih-ey; ++y) {
				for (fx=0; fx<ww; ++fx)
					for (fy=0; fy<wh; ++fy)
						array[fx][fy] = input_image[x+fx-ex][y+fy-ey]; 
				qsort(array[0], ww*wh, sizeof(int), t_compare_ints);
				m_pixel.at(x,y) = array[ww/2][wh/2];
			}

		DeleteImage(array);
	}

	return RRSModule::OK;

};

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new MedianFilter_OMP;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}
