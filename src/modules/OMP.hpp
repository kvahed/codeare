#ifndef _OMP_H_
#define _OMP_H_
#ifdef _OPENMP
#include <omp.h>
#else
inline int  omp_get_thread_num() { return 0;}
inline int  omp_get_num_threads() { return 1;}
inline void omp_set_num_threads(const int) {};
#endif
#endif
