#include "MedianFilter.h"

//integer comparison for qsort
int compare_ints( const void* a, const void* b ) {

	int* arg1 = (int*) a;
	int* arg2 = (int*) b;
	if( *arg1 < *arg2 ) return -1;
	else if( *arg1 == *arg2 ) return 0;
	else return 1;

}

//process the median filter
RRSModule::error_code MedianFilter::ProcessData () {
	
	int m_ww = 25;
	int m_wh = 25;

	int ex = (m_ww / 2), ey = (m_wh / 2);
	int array[m_ww*m_wh];
	
	cout << "MedianFilter::ProcessData() running on " << m_pixel.width() << "x" << 	m_pixel.height() << " image " << endl; 

	int x,y,fx,fy;

	for (x=ex; x<m_pixel.width()-ex; ++x)
		for (y=ey; y<m_pixel.height()-ey; ++y)  {
			for (fx=0; fx<m_ww; ++fx)
				for (fy=0; fy<m_wh; ++fy)
					array[fx+fy*m_ww] = m_pixel.at(x+fx-ex,y+fy-ey);
			qsort(array, m_ww*m_wh, sizeof(int), compare_ints );
			m_pixel.at(x,y) = array[m_ww*m_wh/2] ;
		}

	return RRSModule::OK;

};

