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

#ifndef __NFFT_TRAITS_HPP__
#define __NFFT_TRAITS_HPP__

#include <cmath>
#include <math.h>
#include <stdlib.h>
#include "Complex.hpp"
#include "Vector.hpp"

#include "nfft3.h"

#include "Vector.hpp"

#ifndef USE_NFFT_32_NAMING
#  define nfft_mv_plan_complex mv_plan_complex
#  define nfft_mv_plan_complex mv_plan_complex
#endif

template <class T>
struct NFFTTraits { };

template <> struct NFFTTraits<std::complex<double> > {

    typedef nfft_plan            Plan;    /**< @brief nfft plan (double precision) */
    typedef mri_inh_3d_plan      B0Plan;    /**< @brief nfft plan (double precision) */
    typedef solver_plan_complex  Solver;  /**< @brief nfft solver plan (double precision) */
    typedef double               RT;

    /**
     * @brief            Initialise plan
     *
     * @param  d         Number of dimension (i.e. {1..3})
     * @param  N         Actual dimensions
     * @param  M         Number of k-space samples
     * @param  n         Oversampled N
     * @param  m         Spatial cutoff
     * @param  plan        Forward FT plan
     * @param  solver       Inverse FT plan
     *
     * @return success
     */
    inline static int Init (const Vector<int>& N, size_t M, const Vector<int>& n,
    		int m, Plan& plan, Solver& solver) NOEXCEPT {

        fftw_init_threads();

        Vector<int> _N(N), _n(n);
        int _d (N.size()), _M(M), _m(m);
        unsigned nfft_flags, fftw_flags, solver_flags;

        fftw_flags  = FFTW_ESTIMATE | FFTW_DESTROY_INPUT;
        solver_flags = CGNR | PRECOMPUTE_DAMP | PRECOMPUTE_WEIGHT;
        nfft_flags   = NFFT_OMP_BLOCKWISE_ADJOINT | PRE_PHI_HUT | PRE_PSI |
            MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE;

        nfft_init_guru (&plan, _d, _N.ptr(), _M, _n.ptr(), _m, nfft_flags, fftw_flags);
        solver_init_advanced_complex (&solver, (nfft_mv_plan_complex*) &plan, solver_flags);

        return 0;

    }

    /**
     * @brief            Initialise plan
     *
     * @param  d         Number of dimension (i.e. {1..3})
     * @param  N         Actual dimensions
     * @param  M         Number of k-space samples
     * @param  n         Oversampled N
     * @param  m         Spatial cutoff
     * @param  plan        Forward FT plan
     * @param  solver       Inverse FT plan
     *
     * @return success
     */
    inline static int Init (const Vector<int>& N, size_t M, const Vector<int>& n,
    		int m, RT sigma, B0Plan& plan, Solver& solver) NOEXCEPT {

        Vector<int> _N(N), _n(n);
        int _d (N.size()), _M(M), _m(m);
        unsigned nfft_flags, fftw_flags, solver_flags;

        fftw_flags  = FFTW_ESTIMATE | FFTW_DESTROY_INPUT;
        solver_flags = CGNR | PRECOMPUTE_DAMP | PRECOMPUTE_WEIGHT;
        nfft_flags   = NFFT_OMP_BLOCKWISE_ADJOINT | PRE_PHI_HUT |
                PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE;

        mri_inh_3d_init_guru (&plan, _N.ptr(), _M, _n.ptr(), _m, sigma, nfft_flags, fftw_flags);
        solver_init_advanced_complex (&solver, (nfft_mv_plan_complex*) &plan, solver_flags);

        return 0;

    }



    /**
     * @brief            Inverse FT
     *
     * @param  plan        NFFT plan
     * @param  solver       iNFFT plan
     * @param  maxiter   Maximum NFFT iterations
     * @param  epsilon   Convergence criterium
     *
     * @return           Success
     */
    inline static int ITrafo (Plan& plan, Solver& solver, size_t maxiter = 3,
    		RT epsilon = 3.e-7) NOEXCEPT {

    	/* Initial guess */
    	std::fill_n ((RT*)solver.f_hat_iter, 2*plan.N_total, 0.);

        solver_before_loop_complex(&solver);

        for (size_t l = 0; l < maxiter; l++) {
            if (solver.dot_r_iter < epsilon)
                break;
            solver_loop_one_step_complex(&solver);
        }

        return 0;

    }

    /**
     * @brief            Inverse FT
     *
     * @param  plan        NFFT plan
     * @param  solver       iNFFT plan
     * @param  maxiter   Maximum NFFT iterations
     * @param  epsilon   Convergence criterium
     *
     * @return           Success
     */
    inline static int ITrafo (B0Plan& plan, Solver& solver, size_t maxiter = 3,
    		RT epsilon = 3.e-7) NOEXCEPT {

    	/* Initial guess */
    	std::fill_n ((RT*)solver.f_hat_iter, 2*plan.N_total, 0.);

        solver_before_loop_complex(&solver);

        for (size_t l = 0; l < maxiter; l++) {
            if (solver.dot_r_iter < epsilon)
                break;
            solver_loop_one_step_complex(&solver);
        }

        return 0;

    }


    /**
     * @brief            Forward FT
     *
     * @param  plan        NFFT plan
     *
     * @return           Success
     */
    inline static int Trafo (const Plan& plan) NOEXCEPT {
        nfft_trafo ((Plan*) &plan);
        return 0;
    }

    /**
     * @brief            Forward FT
     *
     * @param  plan        NFFT plan
     *
     * @return           Success
     */
    inline static int Trafo (const B0Plan& plan) NOEXCEPT {
    	mri_inh_3d_trafo ((B0Plan*) &plan);
        return 0;
    }


    /**
     * @brief            Adjoint FT
     *
     * @param  plan        NFFT plan
     *
     * @return           Success
     */
    inline static int Adjoint (const Plan& plan) NOEXCEPT {
        nfft_adjoint ((Plan*) &plan);
        return 0;
    }

    /**
     * @brief            Adjoint FT
     *
     * @param  plan        NFFT plan
     *
     * @return           Success
     */
    inline static int Adjoint  (const B0Plan& plan) NOEXCEPT {
    	mri_inh_3d_adjoint ((B0Plan*) &plan);
        return 0;
    }

    /**
     * @brief            Set weights
     *
     * @param  plan        Plan
     * @param  solver       Solver plan
     * @return           Success
     */
    inline static int  Weights (const Plan& plan, const Solver& solver,
    		const size_t& rank) NOEXCEPT {

    	int err = 0;
        int j, k, z, N = plan.N[0], N2 = N*N, NH = .5*N;
        RT k2, j2, z2;
        
#ifdef NFFT_HAS_MEMBER_FLAGS
        if (solver.flags & PRECOMPUTE_DAMP)
#else
        if (solver.nfft_flags & PRECOMPUTE_DAMP)
#endif
            if (rank == 3) {
                for (j = 0; j < N; ++j) {
                    j2 = j - NH;
                    j2 *= j2;
                    for (k = 0; k < N; ++k) {
                        k2 = k - NH;
                        k2 *= k2;
                        for (z = 0; z < N; ++z) {
                            z2 = z - NH;
                            solver.w_hat[z*N2+j*N+k] = (sqrt(j2+k2+z2*z2) > NH) ? 0. : 1.;
                        }
                    }
                }
            } else {
                for (j = 0; j < N; j++) {
                    j2 = j-NH;
                    j2 *= j2;
                    for (k = 0; k < N; k++) {
                        k2 = k-NH;
                        solver.w_hat[j*N+k] = (sqrt(j2+k2*k2) > NH) ? 0. : 1.;
                    }
                }
            }

        return plan.M_total;

    }


    /**
     * @brief            Precompute PSI
     *
     * @param  plan        NFFT plan
     * @return           Success
     */
    inline static int Psi (Plan& plan) NOEXCEPT {
#ifdef NFFT_HAS_MEMBER_FLAGS
        if(plan.flags & PRE_PSI) /* precompute full psi */
            nfft_precompute_one_psi(&plan);
        if(plan.flags & PRE_FULL_PSI) /* precompute full psi */
            nfft_precompute_full_psi(&plan);
#else
        if(plan.nfft_flags & PRE_PSI) /* precompute full psi */
            nfft_precompute_one_psi(&plan);
        if(plan.nfft_flags & PRE_FULL_PSI) /* precompute full psi */
            nfft_precompute_full_psi(&plan);
#endif
        return 0;
    }



    /**
     * @brief            Finalise plans
     *
     * @param  plan        Plan
     * @param  solver       Solver plan
     * @return           Success
     */
    inline static int Finalize (Plan& plan, Solver& solver) NOEXCEPT {
        solver_finalize_complex(&solver);
        nfft_finalize(&plan);
        return 0;
    }

    /**
     * @brief            Finalise plans
     *
     * @param  plan        Plan
     * @param  solver       Solver plan
     * @return           Success
     */
    inline static int Finalize (B0Plan& plan, Solver& solver) NOEXCEPT {
        solver_finalize_complex(&solver);
        //mri_inh_3d_finalize(&plan);
        return 0;
    }

};

#ifdef HAVE_NFFT3F_UNUSED

template <> struct NFFTTraits<std::complex<float> > {

    typedef nfftf_plan           Plan;    /**< @brief nfft plan (float precision) */
    typedef mrif_inh_3d_plan      B0Plan;    /**< @brief nfft plan (float precision) */
    typedef solverf_plan_complex Solver;  /**< @brief nfft solver plan (float precision) */
    typedef float                RT;

    /**
     * @brief            Initialise plan
     *
     * @param  d         Number of dimension (i.e. {1..3})
     * @param  N         Actual dimensions
     * @param  M         Number of k-space samples
     * @param  n         Oversampled N
     * @param  m         Spatial cutoff
     * @param  plan        Forward FT plan
     * @param  solver       Inverse FT plan
     *
     * @return success
     */
    inline static int Init (const Vector<int>& N, size_t M, const Vector<int>& n,
    		int m, Plan& plan, Solver& solver) NOEXCEPT {

        Vector<int> _N(N), _n(n);
        int _d (N.size()), _M(M), _m(m);
        unsigned nfft_flags, fftw_flags, solver_flags;

        fftw_flags  = FFTW_ESTIMATE | FFTW_DESTROY_INPUT;
        solver_flags = CGNR | PRECOMPUTE_DAMP | PRECOMPUTE_WEIGHT;
        nfft_flags   = NFFT_OMP_BLOCKWISE_ADJOINT | PRE_PHI_HUT |
                PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE;

        nfftf_init_guru (&plan, _d, _N.ptr(), _M, _n.ptr(), _m, nfft_flags, fftw_flags);
        solverf_init_advanced_complex (&solver, (nfftf_mv_plan_complex*) &plan, solver_flags);

        return 0;

    }

    /**
     * @brief            Initialise plan
     *
     * @param  d         Number of dimension (i.e. {1..3})
     * @param  N         Actual dimensions
     * @param  M         Number of k-space samples
     * @param  n         Oversampled N
     * @param  m         Spatial cutoff
     * @param  plan        Forward FT plan
     * @param  solver       Inverse FT plan
     *
     * @return success
     */
    inline static int Init (const Vector<int>& N, size_t M, const Vector<int>& n,
    		int m, RT sigma, B0Plan& plan, Solver& solver) NOEXCEPT {

        Vector<int> _N(N), _n(n);
        int _d (N.size()), _M(M), _m(m);
        unsigned nfft_flags, fftw_flags, solver_flags;

        fftw_flags  = FFTW_ESTIMATE | FFTW_DESTROY_INPUT;
        solver_flags = CGNR | PRECOMPUTE_DAMP | PRECOMPUTE_WEIGHT;
        nfft_flags   = NFFT_OMP_BLOCKWISE_ADJOINT | PRE_PHI_HUT |
                PRE_PSI | MALLOC_X | MALLOC_F_HAT| MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE;

        solverf_init_advanced_complex (&solver, (nfftf_mv_plan_complex*) &plan, solver_flags);

        return 0;

    }



    /**
     * @brief            Inverse FT
     *
     * @param  plan        NFFT plan
     * @param  solver       iNFFT plan
     * @param  maxiter   Maximum NFFT iterations
     * @param  epsilon   Convergence criterium
     *
     * @return           Success
     */
    inline static int ITrafo (Plan& plan, Solver& solver, size_t maxiter = 3,
    		RT epsilon = 3.e-7) NOEXCEPT {

    	/* Initial guess */
    	std::fill_n ((RT*)solver.f_hat_iter, 2*plan.N_total, 0.);

        solverf_before_loop_complex(&solver);

        for (size_t l = 0; l < maxiter; l++) {
            if (solver.dot_r_iter < epsilon)
                break;
            solverf_loop_one_step_complex(&solver);
        }

        return 0;

    }

    /**
     * @brief            Inverse FT
     *
     * @param  plan        NFFT plan
     * @param  solver       iNFFT plan
     * @param  maxiter   Maximum NFFT iterations
     * @param  epsilon   Convergence criterium
     *
     * @return           Success
     */
    inline static int ITrafo (B0Plan& plan, Solver& solver, size_t maxiter = 3,
    		RT epsilon = 3.e-7) NOEXCEPT {

    	/* Initial guess */
    	std::fill_n ((RT*)solver.f_hat_iter, 2*plan.N_total, 0.);

        solverf_before_loop_complex(&solver);

        for (size_t l = 0; l < maxiter; l++) {
            if (solver.dot_r_iter < epsilon)
                break;
            solverf_loop_one_step_complex(&solver);
        }

        return 0;

    }


    /**
     * @brief            Forward FT
     *
     * @param  plan        NFFT plan
     *
     * @return           Success
     */
    inline static int Trafo (const Plan& plan) NOEXCEPT {
        nfftf_trafo ((Plan*) &plan);
        return 0;
    }

    /**
     * @brief            Forward FT
     *
     * @param  plan        NFFT plan
     *
     * @return           Success
     */
    inline static int Trafo (const B0Plan& plan) NOEXCEPT {
    	//mrif_inh_3d_trafo ((B0Plan*) &plan);
        return 0;
    }


    /**
     * @brief            Adjoint FT
     *
     * @param  plan        NFFT plan
     *
     * @return           Success
     */
    inline static int Adjoint (const Plan& plan) NOEXCEPT {
        nfftf_adjoint ((Plan*) &plan);
        return 0;
    }

    /**
     * @brief            Adjoint FT
     *
     * @param  plan        NFFT plan
     *
     * @return           Success
     */
    inline static int Adjoint  (const B0Plan& plan) NOEXCEPT {
    	//mrif_inh_3d_adjoint ((B0Plan*) &plan);
        return 0;
    }

    /**
     * @brief            Set weights
     *
     * @param  plan        Plan
     * @param  solver       Solver plan
     * @return           Success
     */
    inline static int  Weights (const Plan& plan, const Solver& solver,
                                const size_t& rank) NOEXCEPT {

    	int err = 0;
        int j, k, z, N = plan.N[0], N2 = N*N, NH = .5*N;
        RT k2, j2, z2;

#ifdef NFFT_HAS_MEMBER_FLAGS
        if (solver.flags & PRECOMPUTE_DAMP)
#else
        if (solver.nfft_flags & PRECOMPUTE_DAMP)
#endif
            if (rank == 3) {
                for (j = 0; j < N; ++j) {
                    j2 = j - NH;
                    j2 *= j2;
                    for (k = 0; k < N; ++k) {
                        k2 = k - NH;
                        k2 *= k2;
                        for (z = 0; z < N; ++z) {
                            z2 = z - NH;
                            solver.w_hat[z*N2+j*N+k] = (sqrt(j2+k2+z2*z2) > NH) ? 0. : 1.;
                        }
                    }
                }
            } else {
                for (j = 0; j < N; j++) {
                    j2 = j-NH;
                    j2 *= j2;
                    for (k = 0; k < N; k++) {
                        k2 = k-NH;
                        solver.w_hat[j*N+k] = (sqrt(j2+k2*k2) > NH) ? 0. : 1.;
                    }
                }
            }

        return plan.M_total;

    }


    /**
     * @brief            Precompute PSI
     *
     * @param  plan        NFFT plan
     * @return           Success
     */
    inline static int Psi (Plan& plan) NOEXCEPT {
        if(plan.flags & PRE_PSI) /* precompute full psi */
            nfftf_precompute_one_psi(&plan);
        if(plan.flags & PRE_FULL_PSI) /* precompute full psi */
            nfftf_precompute_full_psi(&plan);
        return 0;
    }



    /**
     * @brief            Finalise plans
     *
     * @param  plan        Plan
     * @param  solver       Solver plan
     * @return           Success
     */
    inline static int Finalize (Plan& plan, Solver& solver) NOEXCEPT {
        solverf_finalize_complex(&solver);
        nfftf_finalize(&plan);
        return 0;
    }

    /**
     * @brief            Finalise plans
     *
     * @param  plan        Plan
     * @param  solver       Solver plan
     * @return           Success
     */
    inline static int Finalize (B0Plan& plan, Solver& solver) NOEXCEPT {
        solverf_finalize_complex(&solver);
        //mrif_inh_3d_finalize(&plan);
        return 0;
    }

};

#endif

//
//template<class T> const T Window<T>::KPI = 3.1415926535897932384626433832795028841971693993751;
//template<class T> const T Window<T>::K2PI = 6.2831853071795864769252867665590057683943387987502;
//template<class T> const T Window<T>::KE = 2.7182818284590452353602874713526624977572470937000;
//template<class T> const T Window<T>::P1[] = {
//    1.006897990143384859657820271920512961153421109156614230747188622,
//    0.242805341483041870658834102275462978674549112393424086979586278,
//    0.006898486035482686938510112687043665965094733332210445239567379,
//    0.000081165067173822070066416843139523709162208390998449005642346,
//    4.95896034564955471201271060753697747487292805350402943964e-7,
//    1.769262324717844587819564151110983803173733141412266849e-9,
//    3.936742942676484111899247866083681245613312522754135e-12,
//    5.65030097981781148787580946077568408874044779529e-15,
//    5.267856044117588097078633338366456262960465052e-18,
//    3.111192981528832405775039015470693622536939e-21,
//    1.071238669051606108411504195862449904664e-24,
//    1.66685455020362122704904175079692613e-28
//};
//template<class T> const T Window<T>::Q1[] = {
//    1.000013770640886533569435896302721489503868900260448440877422679934,
//   -0.007438195256024963574139196893944950727405523418354136393367554385,
//    0.000013770655915064256304772604385297068669909609091264440116789601,
//   -1.6794623118559896441239590667288215019925076196457659206142e-8,
//    1.50285363491992136130760477001818578470292828225498818e-11,
//   -1.0383232801211938342796582949062551517465351830706356e-14,
//    5.66233115275307483428203764087829782195312564006e-18,
//   -2.44062252162491829675666639093292109472275754e-21,
//    8.15441695513966815222186223740016719597617e-25,
//   -2.01117218503954384746303760121365911698e-28,
//    3.2919820158429806312377323449729691e-32,
//   -2.70343047912331415988664032397e-36
//};
//template<class T> const T Window<T>::P2[] = {
//    0.4305671332839579065931339658100499864903788418438938270811,
//   -0.2897224581554843285637983312103876003389911968369470222427,
//    0.0299419330186508349765969995362253891383950029259740306077,
//   -0.0010756807437990349677633120240742396555192749710627626584,
//    0.0000116485185631252780743187413946316104574410146692335443,
//   -1.89995137955806752293614125586568854200245376235433e-08
//};
//template<class T> const T Window<T>::Q2[] = {
//    1.0762291019783101702628805159947862543863829764738274558421,
//   -0.7279167074883770739509279847502106137135422309409220238564,
//    0.0762629142282649564822465976300194596092279190843683614797,
//   -0.0028345107938479082322784040228834113914746923069059932628,
//    0.0000338122499547862193660816352332052228449426105409056376,
//   -8.28850093512263912295888947693700479250899073022595e-08
//};
//template<class T> const size_t Window<T>::N1 = sizeof(Window<T>::P1)/sizeof(T);
//template<class T> const size_t Window<T>::M1 = sizeof(Window<T>::Q1)/sizeof(T);
//template<class T> const size_t Window<T>::N2 = sizeof(Window<T>::P2)/sizeof(T);
//template<class T> const size_t Window<T>::M2 = sizeof(Window<T>::Q2)/sizeof(T);
//
//template<class T>
//static inline T chebyshev (const size_t n, const T *c, const T x) NOEXCEPT {
//    T a = c[n-2], b = c[n-1], t;
//    int j;
//    for (size_t j = n - 2; j > 0; j--) {
//        t = c[j-1] - b;
//        b = a + 2.0 * x * b;
//        a = t;
//    }
//    return a + x * b;
//}
//
//template<class T>
//inline static T bessel_i0 (T x) NOEXCEPT {
//    if (x < 0) /* even function */
//        x = -x;
//    if (x == 0.0)
//        return 1.0;
//    if (x <= 15.0) {
//        const T y = x * x;
//        return chebyshev(Window<T>::N1, Window<T>::P1, y)/
//               chebyshev(Window<T>::M1, Window<T>::Q1, y);
//    } else {
//        const T y = (30.0 - x) / x;
//        return (std::exp(x) / std::sqrt(x)) *
//                (chebyshev(Window<T>::N2, Window<T>::P2, y)/
//                 chebyshev(Window<T>::M2, Window<T>::Q2, y));
//    }
//}
//template <class T>
//inline static T phi_hut (T k, int d, const Window<T>& win) NOEXCEPT {
//    return bessel_i0 (win.m * std::sqrt(std::pow(win.b[d], 2.0) -
//            std::pow(2.0 * Window<T>::KPI * k / win.n[d], 2.0)));
//}
//template<class T>
//inline static T phi (T x, int d, const Window<T>& win) NOEXCEPT {
//    return (std::pow(win.m, 2.0) - std::pow(x*win.n[d],2.0) > 0) ?
//            std::sinh(win.b[d] * std::sqrt(std::pow((T)win.m,2.0) -
//                    std::pow(x*win.n[d],2.0)))/(Window<T>::KPI*std::sqrt(std::pow((T)win.m,2.0) -
//                            std::pow(x*win.n[d],2.0))) :
//            (std::pow((T)win.m,2.0) - std::pow((T)x*win.n[d],2.0)) < 0 ?
//                    std::sin(win.b[d] * std::sqrt(std::pow((T)win.n[d]*x,2.0) - std::pow((T)win.m, 2.0))) /
//                    (Window<T>::KPI*std::sqrt(std::pow((T)win.n[d]*x,2.0) - std::pow((T)win.m,2.0))) : 1.0;
//}




#endif //__NFFT_TRAITS_HPP__
