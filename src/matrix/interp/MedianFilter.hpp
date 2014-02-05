#ifndef __MEDIAN_FILTER_HPP__
#define __MEDIAN_FILTER_HPP__

#include "Algos.hpp"
#include "OMP.hpp"
#include "Print.hpp"

inline int
t_compare (const void* a, const void* b) {

    int *ta = (int*) a,
        *tb = (int*) b;
    
	if (*ta < *tb) 
		return -1;
	else if (*ta == *tb) 
		return 0;
	else 
		return 1;

}


template <class T> inline Matrix<T>
medfilt2 (const Matrix<T>& M, const size_t fh = 3, const size_t fw = 3) {

    Matrix<T> ret (vsize(M));

    size_t ih = size(M,0);
    size_t iw = size(M,1);
    
	int x, y, fx, fy, ex = fw/2, ey = fh/2;
    
#pragma omp parallel default (shared) private (y,fx,fy)
    {
        int chunk = iw/omp_get_num_threads();
        Matrix<T> local (fh,fw);
        size_t ne = fw*fh;
        
#pragma omp for schedule(dynamic,chunk)
        for (x=ex; x<iw-ex; ++x)
            for (y=ey; y<ih-ey; ++y) {
                for (fx=0; fx<fw; ++fx)
                    for (fy=0; fy<fh; ++fy)
                        local(fy,fx) = M(y+fy-ey,x+fx-ex); 
                qsort(&local[0], ne, sizeof(T), t_compare);
                ret(y,x) = local[ne/2];
            }
    }
    
    return ret;

}

#endif
