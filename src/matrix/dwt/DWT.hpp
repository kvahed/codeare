/*
 *  codeare Copyright (C) 2013 Daniel Joergens
 *                             Forschungszentrum Juelich, Germany
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
 *
 *  This code follow along the lines of WaveLab 850
 *
 */


# ifndef __DWT_HPP__

# define __DWT_HPP__

/**
 * @brief OMP related makros
 */
# define NUM_THREADS_DWT omp_get_num_threads()
# define OMP_SCHEDULE guided

/**
 * @brief Default wavelet parameters
 */
# define WL_FAM WL_DAUBECHIES
# define WL_MEM 4
# define WL_SCALE 4


/**
 * @brief  Supported wavelet families
 */
enum wlfamily {

	ID = -1,                  /**< Identity transform*/
	WL_DAUBECHIES,
	WL_DAUBECHIES_CENTERED,
	WL_HAAR,
	WL_HAAR_CENTERED,
	WL_BSPLINE,
	WL_BSPLINE_CENTERED

};


# include "Matrix.hpp"
# include "Wavelet.hpp"
# include "Operator.hpp"


/**
 * @brief   Discrete wavelet transform (periodic boundaries) for 2d and 3D case for Matrix template.
 */
template<class T>
class DWT : public Operator<T> {


        /**
         * type definitions
         */
        typedef typename TypeTraits <T> :: RT RT;


    public:


        /**
         * @brief Construct DWT object for images of given side lengths and column major memory scheme.
         *
         * @param  sl1          Side length along first dimension.
         * @param  sl2          Side length along second dimension.
         * @param  sl3          Side length along third dimension.
         * @param  wl_fam       Wavelet family.
         * @param  wl_mem       Member of wavelet family.
         * @param  wl_scale     Decomposition until side length equals 2^wl_scale.
         * @param  num_threads  Number of OMP threads used in parallel regions.
         */
        DWT (const size_t sl1, const size_t sl2, const size_t sl3, const wlfamily wl_fam = WL_FAM,
        		const int wl_mem = WL_MEM, const int wl_scale = WL_SCALE, const int num_threads = NUM_THREADS_DWT) NOEXCEPT
            : _sl1 (sl1),
              _sl2 (sl2),
              _sl3 (sl3),
              _dim (_sl3 == 1 ? 2 : 3),
              _min_sl (_dim == 2 ? std::min (_sl1, _sl2) : std::min (std::min (_sl1, _sl2),_sl3)),
              _min_level (wl_scale),
              _max_level (MaxLevel ()),
              _sl1_scale (_sl1 / pow ((RT)2, (RT)_max_level - _min_level)),
              _sl2_scale (_sl2 / pow ((RT)2, (RT)_max_level - _min_level)),
              _sl3_scale (_sl3 / pow ((RT)2, (RT)_max_level - _min_level)),
              _ld12 (_sl1 * _sl2),
              _wl_fam(wl_fam),
              _fl (wl_mem),
              _modd (_fl/2),
              _meven ((_fl+1)/2),
              _num_threads (num_threads),
              _temp (Vector <T> (_num_threads * std::max (6 * _sl3, std::max (6 * _sl2, 5 * sl1)))),
              dpwt (_dim == 2 ? & DWT <T> :: dpwt2 : & DWT <T> :: dpwt3),
              idpwt (_dim == 2 ? & DWT <T> :: idpwt2 : & DWT <T> :: idpwt3) {
            setupWlFilters <T> (wl_fam, wl_mem, _lpf_d, _lpf_r, _hpf_d, _hpf_r);
        }


        /**
         * @brief       Construct 2D DWT object.
         *
         * @param  sl1          Side length along first dimension.
         * @param  sl2          Side length along second dimension.
         * @param  wl_fam       Wavelet family.
         * @param  wl_mem       Member of wavelet family.
         * @param  wl_scale     Decomposition until side length equals 2^wl_scale.
         * @param  num_threads  Number of OMP threads used in parallel regions.
         */
        DWT (const size_t sl1, const size_t sl2, const wlfamily wl_fam = WL_FAM, const int wl_mem = WL_MEM,
        		const int wl_scale = WL_SCALE, const int num_threads = NUM_THREADS_DWT) NOEXCEPT
        : _sl1 (sl1),
          _sl2 (sl2),
          _sl3 (1),
          _dim (_sl3 == 1 ? 2 : 3),
          _min_sl (_dim == 2 ? std::min (_sl1, _sl2) : std::min (std::min (_sl1, _sl2),_sl3)),
          _min_level (wl_scale),
          _max_level (MaxLevel ()),
          _sl1_scale (_sl1 / pow ((RT)2, (RT)_max_level - _min_level)),
          _sl2_scale (_sl2 / pow ((RT)2, (RT)_max_level - _min_level)),
          _sl3_scale (_sl3 / pow ((RT)2, (RT)_max_level - _min_level)),
          _ld12 (_sl1 * _sl2),
          _wl_fam(wl_fam),
          _fl (wl_mem),
          _modd (_fl/2),
          _meven ((_fl+1)/2),
          _num_threads (num_threads),
          _temp (Vector <T> (_num_threads * std::max (6 * _sl3, std::max (6 * _sl2, 5 * sl1)))),
          dpwt (_dim == 2 ? & DWT <T> :: dpwt2 : & DWT <T> :: dpwt3),
          idpwt (_dim == 2 ? & DWT <T> :: idpwt2 : & DWT <T> :: idpwt3) {
            setupWlFilters <T> (wl_fam, wl_mem, _lpf_d, _lpf_r, _hpf_d, _hpf_r);
        }


        /**
         * @brief       Construct 2D DWT object for square matrices.
         *
         * @param  sl1          Side length along first dimension.
         * @param  wl_fam       Wavelet family.
         * @param  wl_mem       Member of wavelet family.
         * @param  wl_scale     Decomposition until side length equals 2^wl_scale.
         * @param  num_threads  Number of OMP threads used in parallel regions.
         */
        DWT (const size_t sl1,
             const wlfamily wl_fam = WL_FAM, const int wl_mem = WL_MEM, const int wl_scale = WL_SCALE,
             const int num_threads = NUM_THREADS_DWT)
        : _sl1 (sl1),
          _sl2 (_sl1),
          _sl3 (1),
          _dim (_sl3 == 1 ? 2 : 3),
          _min_sl (_dim == 2 ? std::min (_sl1, _sl2) : std::min (std::min (_sl1, _sl2),_sl3)),
          _min_level (wl_scale),
          _max_level (MaxLevel ()),
          _sl1_scale (_sl1 / pow ((RT)2, (RT)_max_level - _min_level)),
          _sl2_scale (_sl2 / pow ((RT)2, (RT)_max_level - _min_level)),
          _sl3_scale (_sl3 / pow ((RT)2, (RT)_max_level - _min_level)),
          _ld12 (_sl1 * _sl2),
          _wl_fam(wl_fam),
          _fl (wl_mem),
          _modd (_fl/2),
          _meven ((_fl+1)/2),
          _num_threads (num_threads),
          _temp (Vector <T> (_num_threads * std::max (6 * _sl3, std::max (6 * _sl2, 5 * sl1)))),
          dpwt (_dim == 2 ? & DWT <T> :: dpwt2 : & DWT <T> :: dpwt3),
          idpwt (_dim == 2 ? & DWT <T> :: idpwt2 : & DWT <T> :: idpwt3) {
            setupWlFilters <T> (wl_fam, wl_mem, _lpf_d, _lpf_r, _hpf_d, _hpf_r);
        }


        virtual
        ~DWT () NOEXCEPT { }


        /**
         * @brief    Forward transform (no constructor calls)
         *
         * @param  m    Signal to decompose
         * @param  res  Resulting DWT
         */
        inline void
        Trafo        (const Matrix <T> & m, Matrix <T> & res) NOEXCEPT {

            assert (   m.Dim (0) == _sl1
                    && m.Dim (1) == _sl2
                    && (_dim == 2 || m.Dim (2) == _sl3)
                    && m.Dim () == res.Dim ());

            /* function pointer */
            (this ->* dpwt) (m, res);

        }


        /**
         * @brief    Adjoint transform (no constructor calls)
         *
         * @param  m    DWT to transform
         * @param  res  Reconstructed signal
         */
        inline void
        Adjoint      (const Matrix <T> & m, Matrix <T> & res) NOEXCEPT {

            assert (   m.Dim (0) == _sl1
                    && m.Dim (1) == _sl2
                    && (_dim == 2 || m.Dim (2) == _sl3)
                    && m.Dim () == res.Dim ());

            /* function pointer */
            (this ->* idpwt) (m, res);

        }


        /**
         * @brief    Forward transform
         *
         * @param  m To transform
         * @return   Transform
         */
        inline Matrix <T>
        operator*    (const Matrix <T> & m) NOEXCEPT {

            if (_wl_fam == ID)
                return m;
            else {
                Matrix <T> res (m);
                Trafo (m, res);
                return res;
            }

        }


        /**
         * @brief    Adjoint transform
         *
         * @param  m To transform
         * @return   Transform
         */
        inline Matrix <T>
        operator->* (const Matrix <T> & m) NOEXCEPT {

            if (_wl_fam == ID)
                return m;
            else {
                Matrix <T> res (m);
                Adjoint (m, res);
                return res;
            }

        }


        inline virtual std::ostream& Print (std::ostream& os) const {
            Operator<T>::Print(os);
            os << "    dims(" << _sl1 << "," << _sl2 << "," << _sl3 << ")" <<
                " WL(" << _wl_fam << "," << _fl << ") threads(" << _num_threads << ")";
            return os;
        }
    

    private:


        /**
         * variable definitions
         */

        // size of valid matrices
        const size_t _sl1;      // side length in first dimension  ('x')
        const size_t _sl2;      // side length in second dimension ('y')
        const size_t _sl3;      // side length in third dimension  ('z')

        // dimension of DWT
        const int _dim;

        const size_t _min_sl;   // minimum side length

        // wavelet scales => (_max_level - _min_level) decompositions
        const int _min_level;   // min. decomposition level
        const int _max_level;   // max. decomposition level

        const size_t _sl1_scale; // starting side length for reconstruction ('x')
        const size_t _sl2_scale; // starting side length for reconstruction ('y')
        const size_t _sl3_scale; // starting side length for reconstruction ('z')
        const int _ld12;        // size of xy - plane

        // wavelet family
        const wlfamily _wl_fam;

        // filter length
        const int _fl;

        // fields for reconstruction in uphi / uplo
        const int _modd;
        const int _meven;

        // number of OMP - threads used in parallel regions
        const int _num_threads;

        // temporary memory used in transform algorithms
        Vector<T> _temp;

        // used transform functions
        void (DWT <T> :: * dpwt) (const Matrix <T> &, Matrix <T> &);
        void (DWT <T> :: * idpwt) (const Matrix <T> &, Matrix <T> &);

        // low pass filters
        RT * _lpf_d;
        RT * _lpf_r;

        // high pass filters
        RT * _hpf_d;
        RT * _hpf_r;


        /**
         * function definitions
         */


        /**
         * @brief           Calculate start level for decomposition.
         *                  (Depends on minimum of side lengths.)
         *
         * @return          Start level.
         */
        inline int
        MaxLevel            () const NOEXCEPT {
            // create vars from mex function
            size_t nn = 1, max_level = 0;
            for (; nn < _min_sl; nn *= 2 )
                max_level ++;
            if (nn  !=  _min_sl) {
                std::cout << "FWT2 requires dyadic length sides" << std::endl;
                assert (false);
            }
            return max_level;
        }


        /**
         *  COPIED FROM WAVELAB IMPLEMENTATION
         */


        /**
         * @brief       Perform forward DWT (periodic boundaries) on 2D data.
         *
         * @param  sig  Signal to be transformed.
         * @param  res  Decomposed signal.
         */
        inline void
        dpwt2		(const Matrix <T> & sig, Matrix <T> & res) NOEXCEPT {

            // assign signal to result matrix
            res = sig;

# pragma omp parallel default (shared), num_threads (_num_threads)
            {

                T * wcplo, * wcphi, * templo, * temphi, * tmp;

                size_t stride;
                int sl1 = _sl1,
                    sl2 = _sl2;
                const int t_num = omp_get_thread_num ();

                // loop over levels of DWT
                for (int j = (_max_level-1); j >= _min_level; --j) {

                    // update stride
                    stride = sl1 * t_num;
                    // update thread's temporary memory address
                    tmp = & _temp [stride];

#pragma omp for //schedule (OMP_SCHEDULE)
                    // loop over lines along first dimension ('columns') of image
                    for (int c2_loc = 0; c2_loc < sl2; c2_loc++) {

                        // access to lowpass part of DWT
                        wcplo = & res [c2_loc * _sl1];
                        // access to highpass part of DWT
                        wcphi = & res [c2_loc * _sl1 + sl1 / 2];

                        // copy part of image to _temp memory
                        copydouble (wcplo, tmp, sl1);

                        // apply low pass filter on current line and write to result matrix
                        downlo (tmp, sl1, wcplo);
                        // apply high pass filter on current line and write to result matrix
                        downhi (tmp, sl1, wcphi);

                    } // loop over lines along first dimension

                    // update stride
                    stride = 2 * sl2 * t_num;
                    // update thread's temporary memory address
                    tmp    = &_temp [stride];
                    templo = &_temp [      sl2 + stride];
                    temphi = &_temp [1.5 * sl2 + stride];

#pragma omp for //schedule (OMP_SCHEDULE)
                    // loop over lines along second dimension ('rows') of image
                    for (int c1_loc = 0; c1_loc < sl1; c1_loc++) {

                        // copy c1-th line of image to temp_mem
                        unpackdouble (& res [0], sl2, _sl1, c1_loc, tmp);

                        // apply low pass filter on current line and write to _temp mem
                        downlo (tmp, sl2, templo);
                        // apply high pass filter on current line and write to _temp mem
                        downhi (tmp, sl2, temphi);

                        // write temp lowpass result to result matrix
                        packdouble (templo, sl2 / 2, _sl1, c1_loc, & res [0]);
                        // write temp highpass result to result matrix
                        packdouble (temphi, sl2 / 2, _sl1, c1_loc, & res [sl2 / 2 * _sl1]);

                    } // loop over lines along second dimension

                    // reduce dimensions for next level
                    sl1 = sl1 / 2;
                    sl2 = sl2 / 2;

                } // loop over levels of DWT

            } // omp parallel

        }


        /**
         * @brief       Perform forward DWT (periodic boundaries) on 3D data.
         *
         * @param  sig  Signal to be transformed.
         * @param  res  Decomposed signal.
         */
        void
        dpwt3       (const Matrix <T> & sig, Matrix <T> & res) NOEXCEPT {

            // assign signal to result matrix
            res = sig;

# pragma omp parallel default (shared), num_threads (_num_threads)
            {

                T * wcplo, * wcphi, * templo, * temphi, * tmp;

                size_t stride;
                int sl1 = _sl1,
                    sl2 = _sl2,
                    sl3 = _sl3;
                const int t_num = omp_get_thread_num ();

                // loop over levels of DWT
                for (int j = (_max_level-1); j >= _min_level; --j) {

                    // update stride
                    stride = sl1 * t_num;
                    // update thread's temporary memory address
                    tmp = & _temp [stride];

# pragma omp for //schedule (OMP_SCHEDULE)
                    // loop over lines along first dimension ('columns') of image
                    for (int c2_loc = 0; c2_loc < sl2 * sl3; c2_loc++) {

                        int c2_glob = (c2_loc / sl2) * _sl1 * _sl2 + (c2_loc % sl2) * _sl1;

                        // access to lowpass part of DWT
                        wcplo = & res [c2_glob /** _sl1*/];
                        // access to highpass part of DWT
                        wcphi = & res [c2_glob /** _sl1*/ + sl1 / 2];

                        // copy part of image to _temp memory
                        copydouble (wcplo, tmp, sl1);

                        // apply low pass filter on current line and write to result matrix
                        downlo (tmp, sl1, wcplo);
                        // apply high pass filter on current line and write to result matrix
                        downhi (tmp, sl1, wcphi);

                    } // loop over lines along first dimension

                    // update stride
                    stride = 2 * sl2 * t_num;
                    // update thread's temporary memory address
                    tmp    = &_temp [stride];
                    templo = &_temp [      sl2 + stride];
                    temphi = &_temp [1.5 * sl2 + stride];

# pragma omp for //schedule (OMP_SCHEDULE)
                    // loop over lines along second dimension ('rows') of image
                    for (int c1_loc = 0; c1_loc < sl1 * sl3; c1_loc++) {

                        int c1_glob = (c1_loc / sl1) * _sl1 * _sl2;

                        // copy c1-th line of image to temp_mem
                        unpackdouble (& res [c1_glob], sl2, _sl1, c1_loc % sl1, tmp);

                        // apply low pass filter on current line and write to _temp mem
                        downlo (tmp, sl2, templo);
                        // apply high pass filter on current line and write to _temp mem
                        downhi (tmp, sl2, temphi);

                        // write temp lowpass result to result matrix
                        packdouble (templo, sl2 / 2, _sl1, c1_loc % sl1, & res [c1_glob]);
                        // write temp highpass result to result matrix
                        packdouble (temphi, sl2 / 2, _sl1, c1_loc % sl1, & res [c1_glob + sl2 / 2 * _sl1]);

                    } // loop over lines along second dimension

                    // update stride
                    stride = 2 * sl3 * t_num;
                    // update thread's temporary memory address
                    tmp    = &_temp [stride];
                    templo = &_temp [      sl3 + stride];
                    temphi = &_temp [1.5 * sl3 + stride];

# pragma omp for //schedule (OMP_SCHEDULE)
                    // loop over lines along third dimension ('third') of image
                    for (int c1_loc = 0; c1_loc < sl1 * sl2; c1_loc++) {

                        int c1_glob = (c1_loc % sl1) + (c1_loc / sl1) * _sl1;

                        // copy c2-th line of image to temp_mem
                        unpackdouble (& res [c1_glob], sl3, _ld12, 0, tmp);

                        // apply low pass filter on current line and write to _temp mem
                        downlo (tmp, sl3, templo);
                        // apply high pass filter on current line and write to _temp mem
                        downhi (tmp, sl3, temphi);

                        // write temp lowpass result to result matrix
                        packdouble (templo, sl3 / 2, _ld12, 0, & res [c1_glob]);
                        // write temp highpass result to result matrix
                        packdouble (temphi, sl3 / 2, _ld12, 0, & res [c1_glob + sl3 / 2 * _ld12]);

                    } // loop over lines along third dimension

                    // reduce dimensions for next level
                    sl1 /= 2;
                    sl2 /= 2;
                    sl3 /= 2;

                } // loop over levels of DWT

            } // omp parallel

        }

        /**
         * @brief           Retrieve scattered data from sequentially stored array.
         *
         * @param  x        Data array to read from.
         * @param  n        Number of data elements to read.
         * @param  stride   Stride between to data elements in context of x.
         * @param  offset   Offset relative to x.
         * @param  y        Data array to be written to.
         */
        inline static void
        unpackdouble	(const T * const x, const int n, const int stride, const int offset, T * const y) NOEXCEPT {
            for (int i = 0; i < n; i++)
                y [i] = x [offset + stride * i];
        }


        /**
         * @brief           Store data using a scattered scheme.
         *
         * @param  x        Data array to read from.
         * @param  n        Number of data elements to store.
         * @param  stride   Stride between to data elements in context of y.
         * @param  offset   Offset relative to y.
         * @param  y        Data array to be written to.
         */
        inline static void
        packdouble		(const T * const x, const int n, const int stride, const int offset, T * const y) NOEXCEPT {
            for (int i = 0; i < n; i++)
                y [offset + stride * i] = x [i];
        }

        /**
         * @brief       Copy data of given length.
         *
         * @param  src  Source.
         * @param  dest Destination.
         * @param  n    Number of data elements.
         */
        inline static void
        copydouble (const T * const src, T * const dest, const int n) NOEXCEPT {
        	std::copy (src, src + n, dest);
        }


        /**
         * @brief       vec(z) = vec(x) + vec(y)
         *
         * @param  x    First summand.
         * @param  y    Second summand.
         * @param  n    Vectors' lengths.
         * @param  z    Result vector.
         */
        inline static void
        adddouble (const T * const x, const T * const y, const int n, T * const z) NOEXCEPT {
            for (int i = 0; i < n; ++i)
                z[i] = x[i] + y[i];
        }


        /**
         * @brief               Perform highpass part of one level forward 1D DWT on signal with given side_length.
         *
         * @param  signal       Signal to be transformed.
         * @param  side_length  Side length of current DWT level.
         * @param  dwt_high     Resulting DWT.
         */
        inline void
        downhi			(const T * const signal, const int side_length, T * const dwt_high) NOEXCEPT {

            int j;
            T s;

            /* highpass version */

            // half of side_length
            const int n2 = side_length / 2;

            // half of filter length
            int mlo = _fl / 2 - 1;

            // adjust lower bound if to low
            if (2 * mlo + 1 - (_fl - 1) < 0)
                mlo++;

            // loop over pixels of dwt_high
            for (int i = mlo; i < n2; i++) {

                // Convolution
                s = 0.;
                for (int h = 0; h < _fl; h++)
                    s += _hpf_d [h] * signal [2 * i + 1 - h];

                dwt_high [i] = s;

            }

            // case: filter_length > side_length => only edge values
            if (mlo > n2)
                mlo = n2;

            /* fix up edge values */

            // loop over edge values
            for (int i = 0; i < mlo; i++) {

                // Convolution
                s = 0.;
                j = 2 * i + 1;

                // loop over filter elements
                for (int h = 0; h < _fl; h++) {

                    // adjust index if it exceeds side_length
                    if (j < 0)
                        j += side_length;
                    s += _hpf_d [h] * signal [j];
                    --j;

                }

                dwt_high [i] = s;

            }
        }


        /**
         * @brief               Perform lowpass part of one level forward 1D DWT on signal with given side_length.
         *
         * @param  signal       Signal to be transformed.
         * @param  side_length  Side length of current DWT level.
         * @param  dwt_low      Resulting DWT.
         */
        void
        downlo  		        (const T * const signal, const int side_length, T * const dwt_low) NOEXCEPT {

            int j;
            T s;

            /*lowpass version */

            // half of side_length (length of dwt_low)
            const int n2 = side_length / 2;
            const int mlo = _fl /2;

            // upper bound for "normal" convolution
            int mhi = n2 - mlo;

            // upper bound to high
            if (2 * mhi + (_fl - 1) >= side_length)
                --mhi;
            // upper bound to low
            if (mhi < 0)
                mhi = -1;

            // loop over pixels of dwt_low
#pragma omp parallel for default (shared) private (s)
            for (int i= 0; i <= mhi; i++) {

                // Convolution
                s = 0.;
                for (int h = 0; h < _fl; h++)
                    s += _lpf_d [h] * signal [2 * i + h];
                dwt_low [i] = s;

            } // loop over pixels of dwt_low


            /* fix up edge values */

            // loop over edge values (periodic boundary)
#pragma omp parallel for default (shared) private (s)
            for (int i = mhi + 1; i < n2; i++) {

                // result of convolution
                s = 0.;

                // start signal index for convolution
                j = 2 * i;

                // loop over filter elements
                for (int h = 0; h < _fl; h++){

                    // adjust index if it exceeds current side_length
                    if (j >= side_length)
                        j -= side_length;
                    s += _lpf_d [h] * signal [j];

                    // update index
                    j++;

                } // perform convolution

                // assign result of convolution
                dwt_low [i] = s;

            }

        }


        /**
         * @brief       Perform inverse DWT (periodic boundaries) on 2D data.
         *
         * @param  wc   Wavelet presentation of 2D data.
         * @param  img  Reconstructed signal.
         */
        inline void
        idpwt2		(const Matrix <T> & wc, Matrix <T> & img) NOEXCEPT {

            // assign dwt to result image
            img = wc;

# pragma omp parallel default (shared) num_threads (_num_threads)
            {

                T * wcplo, * wcphi, * templo, * temphi, * temptop, * tmp;

                size_t stride;
                int sl1 = _sl1_scale,
                    sl2 = _sl2_scale;
                const int t_num = omp_get_thread_num ();

                // loop over levels of backwards DWT
                for (int j = _min_level; j < _max_level; j++) {

                    // update stride
                    stride = 6 * sl2 * t_num;
                    tmp = & _temp [stride];
                    templo  = & _temp [2 * sl2 + stride];
                    temphi  = & _temp [3 * sl2 + stride];
                    temptop = & _temp [4 * sl2 + stride];

# pragma omp for //schedule (OMP_SCHEDULE)
                    // loop over lines along second dimension ('rows') of result image
                    for (int c1_loc = 0; c1_loc < 2 * sl1; c1_loc++)  {

                        // copy lowpass part of current line to temporary memory
                        unpackdouble (& img [0], sl2, _sl1, c1_loc, templo);

                        // copy highpass part of current line to temporary memory
                        unpackdouble (& img [sl2 * _sl1], sl2, _sl1, c1_loc, temphi);

                        // perform lowpass reconstruction
                        uplo (templo, sl2, tmp);
                        // perform highpass reconstruction
                        uphi (temphi, sl2, temptop);

                        // fusion of reconstruction parts
                        adddouble (tmp, temptop, sl2 * 2, tmp);

                        // write back reconstructed line
                        packdouble (tmp, sl2 * 2, _sl1, c1_loc, & img [0]);

                    } // loop over lines along second dimension of result image

                    // update stride
                    stride = 5 * sl1 * t_num;
                    tmp = & _temp [stride];
                    templo = & _temp [    sl1 + stride];
                    temphi = & _temp [3 * sl1 + stride];

# pragma omp for //schedule (OMP_SCHEDULE)
                    // loop  over lines along first dimension ('columns') of result image
                    for (int c2_loc = 0; c2_loc < 2 * sl2; c2_loc++) {

                        // assign address of current line's lowpass part
                        wcplo = & img [c2_loc * _sl1];
                        // assign address of current line's highpass part
                        wcphi = & img [c2_loc * _sl1 + sl1];

                        // copy lowpass part to temporary memory
                        copydouble (wcplo, tmp, sl1);

                        // perform lowpass reconstruction
                        uplo (wcplo, sl1, templo);
                        // perform highpass reconstruction
                        uphi (wcphi, sl1, temphi);

                        // combine reconstructed parts and write back to current line
                        adddouble (templo, temphi, sl1 * 2, wcplo);

                    } // loop over lines along first dimension ('columns') of result image

                    // update current row / column size
                    sl2 *= 2;
                    sl1 *= 2;

                } // loop over levels of backwards DWT

            } // omp parallel

        }


        /**
         * @brief       Perform inverse DWT (periodic boundaries) on 3D data.
         *
         * @param  wc   Wavelet presentation of 3D data.
         * @param  img  Reconstructed signal.
         */
        inline  void idpwt3      (const Matrix <T> & wc, Matrix <T> & img)  NOEXCEPT {

            // assign dwt to result image
            img = wc;

# pragma omp parallel default (shared) num_threads (_num_threads)
            {

                T * wcplo, * wcphi, * templo, * temphi, * temptop, * tmp;

                size_t stride;
                int sl1 = _sl1_scale,
                    sl2 = _sl2_scale,
                    sl3 = _sl3_scale;
                const int t_num = omp_get_thread_num ();

                // loop over levels of backwards DWT
                for (int j = _min_level; j < _max_level; j++) {

                    // update stride
                    stride = 6 * sl3 * t_num;
                    tmp = & _temp [stride];
                    templo  = & _temp [2 * sl3 + stride];
                    temphi  = & _temp [3 * sl3 + stride];
                    temptop = & _temp [4 * sl3 + stride];

# pragma omp for //schedule (OMP_SCHEDULE)
                    // loop over lines along third dimension ('third') of result image
                    for (int c1_loc = 0; c1_loc < 2 * sl1 * 2 * sl2; c1_loc++)
                    {

                        int c1_glob = (c1_loc % (2 * sl1)) + (c1_loc / (2 * sl1)) * _sl1;

                        // copy lowpass part of current line to temporary memory
                        unpackdouble (& img [c1_glob], sl3, _ld12, 0, templo);

                        // copy highpass part of current line to temporary memory
                        unpackdouble (& img [c1_glob + sl3 * _ld12], sl3, _ld12, 0, temphi);

                        // perform lowpass reconstruction
                        uplo (templo, sl3, tmp);
                        // perform highpass reconstruction
                        uphi (temphi, sl3, temptop);

                        // fusion of reconstruction parts
                        adddouble (tmp, temptop, sl3 * 2, tmp);

                        // write back reconstructed line
                        packdouble (tmp, sl3 * 2, _ld12, 0, & img [c1_glob]);

                    } // loop over lines along third dimension of result image

                    // update stride
                    stride = 6 * sl2 * t_num;
                    tmp = & _temp [stride];
                    templo  = & _temp [2 * sl2 + stride];
                    temphi  = & _temp [3 * sl2 + stride];
                    temptop = & _temp [4 * sl2 + stride];

# pragma omp for //schedule (OMP_SCHEDULE)
                    // loop over lines along second dimension ('rows') of result image
                    for (int c1_loc = 0; c1_loc < 2 * sl1 * 2 * sl3; c1_loc++)
                    {

                        int c1_glob = (c1_loc / (2 * sl1)) * _sl1 * _sl2;

                        // copy lowpass part of current line to temporary memory
                        unpackdouble (& img [c1_glob], sl2, _sl1, c1_loc % (2 * sl1), templo);

                        // copy highpass part of current line to temporary memory
                        unpackdouble (& img [c1_glob + sl2 * _sl1], sl2, _sl1, c1_loc % (2 * sl1), temphi);

                        // perform lowpass reconstruction
                        uplo (templo, sl2, tmp);
                        // perform highpass reconstruction
                        uphi (temphi, sl2, temptop);

                        // fusion of reconstruction parts
                        adddouble (tmp, temptop, sl2 * 2, tmp);

                        // write back reconstructed line
                        packdouble (tmp, sl2 * 2, _sl1, c1_loc % (2 * sl1), & img [c1_glob]);

                    } // loop over lines along second dimension of result image

                    // update stride
                    stride = 5 * sl1 * t_num;
                    tmp = & _temp [stride];
                    templo = & _temp [    sl1 + stride];
                    temphi = & _temp [3 * sl1 + stride];

# pragma omp for //schedule (OMP_SCHEDULE)
                    // loop  over lines along first dimension ('columns') of result image
                    for (int c2_loc = 0; c2_loc < 2 * sl2 * 2 * sl3; c2_loc++)
                    {

                        int c2_glob = (c2_loc / (2 * sl2)) * _sl2 * _sl1 + (c2_loc % (2 * sl2)) * _sl1;

                        // assign address of current line's lowpass part
                        wcplo = & img [c2_glob];
                        // assign address of current line's highpass part
                        wcphi = & img [c2_glob + sl1];

                        // copy lowpass part to temporary memory
                        copydouble (wcplo, tmp, sl1);

                        // perform lowpass reconstruction
                        uplo (wcplo, sl1, templo);
                        // perform highpass reconstruction
                        uphi (wcphi, sl1, temphi);

                        // combine reconstructed parts and write back to current line
                        adddouble (templo, temphi, sl1 * 2, wcplo);

                    } // loop over lines along first dimension ('columns') of result image

                    // update current row / column size
                    sl2 *= 2;
                    sl1 *= 2;
                    sl3 *= 2;

                } // loop over levels of backwards DWT

            } // omp parallel

        }


        /**
         * @brief               Perform lowpass reconstruction.
         *
         * @param  wc           1D wavelet representation.
         * @param  side_length  Length of wc.
         * @param  signal       Reconstructed signal.
         */
        inline void
        uplo		(const T * const wc, const int side_length, T * const signal) NOEXCEPT {

            int j;
            T s, s_odd;

            /*lowpass version */

            /* away from edges */

            // loop over regular signal indices
            for (int i = _meven; i < side_length; i++)
            {

                // init convolution results
                s = 0.;
                s_odd = 0.;

                // perform convolution for even and odd filter indices
                for (int h = 0; h < _modd; h++)
                {

                    // even filter index
                    s += _lpf_r [2 * h] * wc [i - h];
                    // odd filter index
                    s_odd += _lpf_r [2 * h + 1] * wc [i - h];

                }
                // case of odd filter_length (-> more even indices: start with index 0)
                if (_meven > _modd)
                    s += _lpf_r [2 * _meven] * wc [i - _meven];

                // assign convolution results
                signal [2 * i] = s;
                signal [2 * i + 1] = s_odd;

            } // loop over regular signal indices


            /* fix up edge values */

            // upper bound for filter indices
            int mmax = _meven;
            // possible correction if mmax greater than current side length
            if (mmax > side_length)
                mmax = side_length;

            // loop over edge values
            for (int i = 0; i < mmax; i++)
            {

                // init convolution results
                s = 0.;
                s_odd = 0.;
                // set start index of wavelet coefficients
                j = i;

                 // perform convolution
                for (int h = 0; h < _modd; h++)
                {

                    // correct current wavelet coeff's index if needed
                    if (j < 0)
                        j += side_length;

                    // even part of convolution
                    s += _lpf_r [2 * h] * wc [j];
                    // odd part of convolution
                    s_odd += _lpf_r [2 * h + 1] * wc [j];

                    // update index
                    --j;

                } // perform convolution

                // case of odd filter_length
                if (_meven > _modd)
                    s += _lpf_r [2 * _meven] * wc [j];

                // assign convolution results
                signal [2 * i] = s;
                signal [2 * i + 1] = s_odd;

            } // loop over edge values

        }


        /**
         * @brief               Perform highpass reconstruction.
         *
         * @param  wc           1D wavelet representation.
         * @param  side_length  Length of wc.
         * @param  signal       Reconstructed signal.
         */
        inline void
        uphi		(const T * const wc, const int side_length, T * const signal) NOEXCEPT {

            int j;
            T s, s_odd;

            /*hipass version */

            /* away from edges */

            // loop over regular signal indices
            for (int i = 0; i < side_length - _meven; i++){

                // init convolution results
                s = 0.;
                s_odd = 0.;

                // perform convolution for even and odd filter indices
                for (int h = 0; h < _modd; h++) {

                    // even filter index
                    s += _hpf_r [2 * h] * wc [i + h];
                    // odd filter index
                    s_odd += _hpf_r [2 * h + 1] * wc [i + h];

                } // perform convolution

                // case of odd filter_length
                if (_meven > _modd)
                    s += _hpf_r [2 * _meven] * wc [i + _meven];

                // assign convolution results
                signal [2 * i + 1] = s;
                signal [2 * i] = s_odd;

            } // loop over regular signal indices


            /* fix up edge values */

            // lower bound for indices of edge values
            int mmin = side_length - _meven;
            // possible correction if mmin less than zero
            if (mmin < 0)
                mmin = 0;

            // loop over edge values
            for (int i = mmin; i < side_length; i++) {

                // init convolution results
                s = 0.;
                s_odd = 0.;
                // start index of wavelet coefficients
                j = i;

                // perform convolution for even and odd indices
                for (int h = 0; h < _meven; h++) {

                    // correct current wavelet coeff's index if needed
                    if (j >= side_length)
                        j -= side_length;

                    // even filter index
                    s += _hpf_r [2 * h] * wc [j];
                    // odd filter index
                    s_odd += _hpf_r [2 * h + 1] * wc [j];

                    // update current wc index
                    j++;

                } // perform convolution

                // assign convolution results
                signal [2 * i + 1] = s;
                signal [2 * i] = s_odd;

            } // loop over edge values

        }

};


# endif // __DWT_HPP__
