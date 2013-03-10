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

#include "Workspace.hpp"
#include "Params.hpp"
#include "OMP.hpp"

#include <fftw3.h>

template <class T>
struct FTTraits { };

/**
 * @brief C++ friendly interface to complex FFTW (single precision)
 */
template<>
struct FTTraits<float> {
	
	typedef fftwf_plan    Plan; /**< @brief fftw plan (float precision) */
	typedef fftwf_complex Type; /**< @brief fftw complex data type (float precision) */
	

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
	DFTPlan (int rank, const int *n, Type *in, Type *out, int sign, unsigned flags) {
		InitThreads();
	    return fftwf_plan_dft (rank, n, in, out, sign, flags);
	}
	

		/**
	 * @brief         Initialise FFTWF threads (Should only be done once)
	 *
	 * @return        Success
	 */
	static bool 
	InitThreads () {

		Workspace& ws = Workspace::Instance();
		int nt = boost::any_cast<int>(ws.p["FFTWFThreads"]);

		if (nt > 0)
			return true;

#pragma omp parallel
		{
			nt = omp_get_num_threads();
			omp_set_dynamic(false);
		}		

		printf ("  Initialising %d-threaded FFTWF ... ", nt);
		fflush (stdout);

#ifdef HAVE_FFTWF_THREADS
		bool ok = fftwf_init_threads();
		if (ok) {
			PlanWithNThreads (nt);
			ws.PSet("FFTWFThreads", nt);
		}
#endif

		printf ("done\n");

		return true;
		
	}
	

	/**
	 * @brief        Create plans with N threads
	 * 
	 * @param  nt    # of threads
	 */
	static void 
	PlanWithNThreads (const int& nt) { 
		fftwf_plan_with_nthreads (nt); 
	}


	/**
	 * @brief        Inlined memory allocation for performance
	 *
	 * @param  n     # of elements
	 * @return       Memory address
	 */
	static inline Type*
	Malloc (size_t n) {
		return (Type*) fftwf_malloc (n * sizeof(Type));
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
	
};


/**
 * @brief C++ friendly interface to complex FFTW (double precision)
 */
template<>
struct FTTraits<double> {
	
	typedef fftw_plan    Plan;  /**< @brief fftw plan (double precision) */
	typedef fftw_complex Type;  /**< @brief fftw complex data type (double precision) */
	

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
	DFTPlan (int rank, const int *n, Type *in, Type *out, int sign, unsigned flags) {                    
	    return fftw_plan_dft (rank, n, in, out, sign, flags);
	}
	

	/**
	 * @brief         Initialise FFTW threads (Should only be done once)
	 *
	 * @return        Success
	 */
	static bool 
	InitThreads () {

		Workspace& ws = Workspace::Instance();
		int nt = boost::any_cast<int>(ws.p["FFTWFThreads"]);

		if (nt > 0)
			return true;

#pragma omp parallel
		{
			nt = omp_get_num_threads();
		}		

		printf ("  Initialising %d-threaded FFTW ... ", nt);
		fflush (stdout);

#ifdef HAVE_FFTWF_THREADS
		bool ok = fftw_init_threads();
		if (ok) {
			PlanWithNThreads (nt);
			ws.PSet("FFTWThreads", nt);
		}
#endif

		printf ("done\n");

		return true;
		
	}
	

	/**
	 * @brief        Create plans with N threads
	 * 
	 * @param   nt   # of threads
	 */
	static void 
	PlanWithNThreads (const int& nt) { 
		fftw_plan_with_nthreads (nt); 
	}


	/**
	 * @brief        Inlined memory allocation for performance
	 *
	 * @param  n     # of elements
	 * @return       Memory address
	 */
	static inline Type*
	Malloc (size_t n) {
		return (Type*) fftw_malloc (n * sizeof(Type));
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

};

#endif
