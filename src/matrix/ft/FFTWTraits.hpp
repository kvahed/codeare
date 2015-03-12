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

#ifndef __FFTW_TRAITS_HPP__
#define __FFTW_TRAITS_HPP__

//#include "Workspace.hpp"
#include "Params.hpp"
#include "OMP.hpp"

#include <fftw3.h>

static Params p;

template <class T>
struct FTTraits { };

/**
 * @brief C++ friendly interface to complex FFTW (single precision)
 */
template<>
struct FTTraits<float> {
	
	typedef fftwf_plan    Plan; /**< @brief fftw plan (float precision) */
	typedef fftwf_complex T; /**< @brief fftw complex data type (float precision) */
	typedef float         otype;
	

	/**
	 * @brief         DFT plan
	 *
	 * @param  rank   FT dimesionality
	 * @param  n      Size lengths of individual dimensions
	 * @param  in     Input memory
	 * @param  out    Output memory
	 * @param  sign   FT direction
	 * @param  flags  FFTW flags
	 *
	 * @return        Plan
	 */
	static inline Plan 
	DFTPlan (int rank, const int *n, T *in, T *out, int sign, unsigned flags) {
		InitThreads();
#ifdef _OPENMP
		fftwf_import_wisdom_from_filename("codeare_threads.plan");
#else
		fftwf_import_wisdom_from_filename("codeare_single.plan");
#endif
		Plan plan = fftwf_plan_dft (rank, n, in, out, sign, flags);
#ifdef _OPENMP
		fftwf_export_wisdom_to_filename("codeare_threads.plan");
#else
		fftwf_export_wisdom_to_filename("codeare_single.plan");
#endif
		return plan;
	}
	

	/**
	 * @brief         Initialise FFTWF threads (Should only be done once)
	 *
	 * @return        Success
	 */
	static inline bool
	InitThreads () {

		if (p.exists("FFTWThreads"))
			return true;

		int nt;

#pragma omp parallel
		{
			nt = omp_get_num_threads();
		}		

#ifdef _OPENMP
		if (fftwf_init_threads()) {
			fftwf_plan_with_nthreads (nt);
			p["FFTWThreads"] = nt;
		}
#endif

		return true;
		
	}

	
	static inline Plan DFTPlanMany (int rank, const int* n, int howmany, T* in, T* out, int dir) {
		InitThreads();
		return fftwf_plan_many_dft (rank, n, howmany, in, NULL, 1, *n, out, NULL, 1, *n, dir, FFTW_ESTIMATE);
	}


	/**
	 * @brief        Inlined memory allocation for performance
	 *
	 * @param  n     # of elements
	 * @return       Memory address
	 */
	static inline T*
	Malloc (size_t n) {
		return (T*) fftwf_malloc (n * sizeof(T));
	}
	

	/**
	 * @brief        Free allocated memory
	 *
	 * @param  p     Memory address
	 */
	static inline void 
	Free (void* p) {
		fftwf_free (p);
	}
	

	/**
	 * @brief        Clean up tools
	 */
	static inline void 
	CleanUp () {
#ifdef _OPENMP
		fftwf_cleanup_threads();
#endif
		fftwf_cleanup ();
	}
	

	/**
	 * @brief        Destroy plan
	 *
	 * @param  p     Plan to be destroyed
	 */
	static inline void 
	Destroy (Plan p) {
		fftwf_destroy_plan (p); 
	}
	

	/**
	 * @brief        Execute plan
	 *
	 * @param  p     Plan to be executed
	 */
	static inline void 
	Execute (Plan p) {
		fftwf_execute (p); 
	}        
	
	/**
	 * @brief        Execute plan
	 *
	 * @param  p     Plan to be executed
	 */
	static inline void
	Execute (Plan p, T* in, T* out) {
		fftwf_execute_dft (p, in, out);
	}

};


/**
 * @brief C++ friendly interface to complex FFTW (double precision)
 */
template<>
struct FTTraits<double> {
	
	typedef fftw_plan    Plan;  /**< @brief fftw plan (double precision) */
	typedef fftw_complex T;  /**< @brief fftw complex data type (double precision) */
	typedef double      otype;

	/**
	 * @brief         DFT plan
	 *
	 * @param  rank   FT dimesionality
	 * @param  n      Size lengths of individual dimensions
	 * @param  in     Input memory
	 * @param  out    Output memory
	 * @param  sign   FT direction
	 * @param  flags  FFTW flags
	 *
	 * @return        Plan
	 */
	static inline Plan 
	DFTPlan (int rank, const int *n, T *in, T *out, int sign, unsigned flags) {
		InitThreads();
#ifdef _OPENMP
		fftw_import_wisdom_from_filename("codeare_threads.plan");
#else
		fftw_import_wisdom_from_filename("codeare_single.plan");
#endif
		Plan plan = fftw_plan_dft (rank, n, in, out, sign, flags);
#ifdef _OPENMP
		fftw_export_wisdom_to_filename("codeare_threads.plan");
#else
		fftw_export_wisdom_to_filename("codeare_single.plan");
#endif
		return plan;
	}
	

	/**
	 * @brief         Initialise FFTW threads (Should only be done once)
	 *
	 * @return        Success
	 */
	static inline bool
	InitThreads () {

		if (p.exists("FFTWThreads"))
			return true;

		int nt;

#pragma omp parallel
		{
			nt = omp_get_num_threads();
		}		

#ifdef _OPENMP
		if (fftw_init_threads()) {
			fftw_plan_with_nthreads (nt);
			p["FFTWThreads"] = nt;
		}
#endif

		return true;
		
	}
	

	static inline Plan DFTPlanMany (int rank, const int* n, int howmany, T* in, T* out, int dir) {
		InitThreads();
		return fftw_plan_many_dft (rank, n, howmany, in, NULL, 1, *n, out, NULL, 1, *n, dir, FFTW_ESTIMATE);
	}

	/**
	 * @brief        Inlined memory allocation for performance
	 *
	 * @param  n     # of elements
	 * @return       Memory address
	 */
	static inline T*
	Malloc (size_t n) {
		return (T*) fftw_malloc (n * sizeof(T));
	}
	

	/**
	 * @brief        Free allocated memory
	 *
	 * @param  p     Memory address
	 */
	static inline void
	Free (void* p) {
		fftw_free(p);
	}


	/**
	 * @brief        Clean up tools
	 */
	static inline void 
	CleanUp () {
#ifdef _OPENMP
		fftw_cleanup_threads();
#endif
		fftw_cleanup ();
	}
	

	/**
	 * @brief        Destroy plan
	 *
	 * @param   p    Plan to be destroyed
	 */
	static inline void 
	Destroy (Plan p) {
		fftw_destroy_plan (p); 
	}


	/**
	 * @brief       Execute plan
	 *
	 * @param  p    Plan to be executed
	 */
	static inline void
	Execute (Plan p) {
		fftw_execute (p); 
	}        

	/**
	 * @brief        Execute plan
	 *
	 * @param  p     Plan to be executed
	 */
	static inline void
	Execute (Plan p, T* in, T* out) {
		fftw_execute_dft (p, in, out);
	}

};

#endif
