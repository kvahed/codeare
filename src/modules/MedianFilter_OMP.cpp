/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
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

#include "MedianFilter_OMP.hpp"

using namespace RRStrategy;

template <class T> T** CreateImage (int n, int m) {

	T** pp = new T* [n];
	pp[0]  = new T [n*m];

	for (int i = 1; i < n; ++i) 
		pp[i] = pp[i-1] + m;

	return pp;

}

int t_compare_ints 
(const void* a, const void* b) {

	int* arg1 = (int*) a;
	int* arg2 = (int*) b;

	if ( *arg1 < *arg2 ) 
		return -1;
	else if( *arg1 == *arg2 ) 
		return 0;
	else 
		return 1;

}

template <class T> void 
DeleteImage (T** pp) {

	delete [] pp[0];
	delete [] pp;

}


error_code
MedianFilter_OMP::Init () {

	printf ("Intialising MedianFilter ...\n");

    double temp;
    
    Attribute ("wh", &temp);
    m_wh = (unsigned short)temp;
    printf ("  Filter size : %ix", m_wh);

    Attribute ("ww", &temp);
    m_ww = (unsigned short)temp;
    printf ("%i\n", m_ww);
    
    m_uname = Attribute ("uname");
    printf ("  Unique name : %s\n", m_uname);    

	printf ("... done.\n");


}

error_code 
MedianFilter_OMP::Process () {
	
	ticks cgstart = getticks();
	printf ("Processing MedianFilter ...\n");

	Matrix<short>& img = Get<short> (m_uname);
    if (img.Size() == 1)
        return OK;

	int iw = img.Width(); //image width 
	int ih = img.Height(); //image height
	
	int ex = (m_ww / 2), ey = (m_wh / 2);
	
	int x,y,fx,fy;
	int** input_image  = CreateImage<int>(iw,ih);
	
#pragma omp parallel default(shared)  private(y,fx,fy)
	{
		int tid      = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int chunk    = iw/nthreads; //chunk-size for loop splitting
		
		if (tid==0) 
			std::cout << "MedianFilter_OMP::Process() running on " << iw << "x" << ih 
					  << " image with " << nthreads << " threads" << std::endl; 

#pragma omp for schedule(dynamic,chunk)
		for (x=0; x<iw; ++x)
			for (y=0; y<ih; ++y)
				input_image[x][y]  = img(x,y);

		int** array = CreateImage<int>(m_ww,m_wh); //local to each thread !

#pragma omp for  schedule(dynamic,chunk)
		for (x=ex; x<iw-ex; ++x)
			for (y=ey; y<ih-ey; ++y) {
				for (fx=0; fx<m_ww; ++fx)
					for (fy=0; fy<m_wh; ++fy)
						array[fx][fy] = input_image[x+fx-ex][y+fy-ey]; 
				qsort(array[0], m_ww*m_wh, sizeof(int), t_compare_ints);
				img(x,y) = array[m_ww/2][m_wh/2];
			}

		DeleteImage(array);
	}

	printf ("... done. WTime: %.4f seconds.\n\n", elapsed(getticks(), cgstart) / Toolbox::Instance()->ClockRate());

	return OK;

}

// the class factories
extern "C" DLLEXPORT ReconStrategy* create  ()                 {
    return new MedianFilter_OMP;
}

extern "C" DLLEXPORT void           destroy (ReconStrategy* p) {
    delete p;
}
