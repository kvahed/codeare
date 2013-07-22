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
 * OMP related makros
 */
# define NUM_THREADS_DWT 4
# define OMP_SCHEDULE guided


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



/**
 * @brief 2D Discrete wavelet transform for Matrix template (from GSL)
 */
template<class T>
class DWT {


    public:


        /**
         * @brief Construct 2D Wavelet transform with wavelet class and side length
         *
         * @param  sl      Side length
         */
        DWT (const size_t sl1, const size_t sl2, const wlfamily wl_fam = WL_DAUBECHIES, const int wl_mem = 4,
        		const int wl_scale = 4, const int num_threads = NUM_THREADS_DWT, const int dim = 2)
            : _dim (dim),
              _sl1 (sl1),
              _sl2 (sl2),
              _sl3 (1),
              _fl (wl_mem),
              temp (container<T>(num_threads * MAX (6 * sl2, 5 * sl1))),
              _wl_scale (wl_scale),
              _fam(wl_fam),
              _J (0),
              _num_threads (num_threads) {

            setupWlFilters <T> (wl_fam, wl_mem, _lpf_d, _lpf_r, _hpf_d, _hpf_r);


            if (_dim == 2)
                _min_sl = MIN (_sl1, _sl2);
            else
                _min_sl = MIN (MIN (_sl1, _sl2),_sl3);

            // create vars from mex function
            int nn;
            for (nn = 1; nn < _min_sl; nn *= 2 )
                _J ++;
            if (nn  !=  _min_sl){
                std::cout << "FWT2 requires dyadic length sides" << std::endl;
                assert (false);
            }

            // calc 2^...
            int tmp_two = 1;
            for (nn = _wl_scale; nn < _J; nn++)
                tmp_two *= 2;

//            std::cout << " min_sl: " << _min_sl << std::endl;
            _sl1_scale = _sl1 / tmp_two;
            _sl2_scale = _sl2 / tmp_two;
            _sl3_scale = _sl3 / tmp_two;

        }

        virtual
        ~DWT ()
        { }


        /**
         * @brief    Forward transform (no constructor calls)
         *
         * @param  m    Signal to decompose
         * @param  res  Resulting DWT
         */
        inline
        void
        Trafo        (const Matrix <T> & m, Matrix <T> & res)
        {

            assert (   m.Dim (0) == _sl1
                    && m.Dim (1) == _sl2
                    && (_dim == 2 || m.Dim (2) == _sl3)
                    && m.Dim () == res.Dim ());

            // call dpwt2
            dpwt2 (m, _wl_scale, _J, temp, res);

        }


        /**
         * @brief    Adjoint transform (no constructor calls)
         *
         * @param  m    DWT to transform
         * @param  res  Reconstructed signal
         */
        inline
        void
        Adjoint      (const Matrix <T> & m, Matrix <T> & res)
        {

            assert (   m.Dim (0) == _sl1
                    && m.Dim (1) == _sl2
                    && (_dim == 2 || m.Dim (2) == _sl3)
                    && m.Dim () == res.Dim ());

            // call idpwt2
            idpwt2 (m, _wl_scale, _J, temp, res);

        }


        /**
         * @brief    Forward transform
         *
         * @param  m To transform
         * @return   Transform
         */
        inline
        Matrix <T>
        operator*    (const Matrix <T> & m) {

            if (_fam == ID)
                return m;
            else
            {
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
        inline
        Matrix <T>
        operator->* (const Matrix <T> & m) {

            if (_fam == ID)
                return m;
            else
            {
                Matrix <T> res (m);
                Adjoint (m, res);
                return res;
            }

        }


    private:


        /**
         * type definitions
         */
        typedef typename TypeTraits <T> :: RT RT;


        /**
         * variable definitions
         */

        // dimension of DWT
        const int _dim;

        // low pass filters
        RT * _lpf_d;
        RT * _lpf_r;

        // high pass filters
        RT * _hpf_d;
        RT * _hpf_r;

        // maximum size of matrix
        const size_t _sl1; // side length in first dimension  ('x')
        const size_t _sl2; // side length in second dimension ('y')
        const size_t _sl3; // side length in third dimension  ('z')
        size_t _min_sl;
        size_t _sl1_scale;
        size_t _sl2_scale;
        size_t _sl3_scale;

        int _J;

        // filter length
        const size_t _fl;

        // temporary memory
        container<T> temp;

        // wavelet scale (max. decomposition level)
        const int _wl_scale;

        const wlfamily _fam;

        const int _num_threads;


        /**
         * function definitions
         */


        /**
         *  COPIED FROM WAVELAB IMPLEMENTATION
         */


        void
        dpwt2		(const Matrix <T> & sig, const int ell, const int J,
             		 container <T> & temp_mem, Matrix <T> & res)
        {

            T * wcplo, * wcphi, * templo, * temphi, * tmp;

            // assign signal to result matrix
            res = sig;

# pragma omp parallel default (shared), private (wcplo, wcphi, temphi, templo, tmp)\
                                        num_threads (_num_threads)
            {

                size_t stride;
                int sl1 = _sl1,
                    sl2 = _sl2,
                    sl3 = _sl3;
                const int t_num = omp_get_thread_num ();
                int c1_glob,
                    c2_glob,
                    c3_glob;

                // loop over levels of DWT
                for (int j = (J-1); j >= ell; --j)
                {

                    // update stride
                    stride = sl1 * t_num;
                    // update thread's temporary memory address
                    tmp = & temp_mem [stride];

# pragma omp for schedule (OMP_SCHEDULE)
                    // loop over lines along first dimension ('columns') of image
                    for (int c2_loc = 0; c2_loc < sl2 * sl3; c2_loc++)
                    {

//                        c2_glob = c2 % _sl2 +

                        // access to lowpass part of DWT
                        wcplo = & res [c2_loc * _sl1];
                        // access to highpass part of DWT
                        wcphi = & res [c2_loc * _sl1 + sl1 / 2];

                        // copy part of image to temp memory
                        copydouble (wcplo, tmp, sl1);

                        // apply low pass filter on current line and write to result matrix
                        downlo (tmp, sl1, wcplo);
                        // apply high pass filter on current line and write to result matrix
                        downhi (tmp, sl1, wcphi);

                    } // loop over lines along first dimension

                    // update stride
                    stride = 2 * sl2 * t_num;
                    // update thread's temporary memory address
                    tmp = & temp_mem [stride];
                    templo = & temp_mem [      sl2 + stride];
                    temphi = & temp_mem [1.5 * sl2 + stride];

# pragma omp for schedule (OMP_SCHEDULE)
                    // loop over lines along second dimension ('rows') of image
                    for (int c1 = 0; c1 < sl1 * sl3; c1++)
                    {

                        // copy c1-th line of image to temp_mem
                        unpackdouble (res.Memory (0), sl2, _sl1, c1, tmp);

                        // apply low pass filter on current line and write to temp mem
                        downlo (tmp, sl2, templo);
                        // apply high pass filter on current line and write to temp mem
                        downhi (tmp, sl2, temphi);

                        // write temp lowpass result to result matrix
                        packdouble (templo, sl2 / 2, _sl1, c1, & res [0]);
                        // write temp highpass result to result matrix
                        packdouble (temphi, sl2 / 2, _sl1, c1, & res [sl2 / 2 * _sl1]);

                    } // loop over lines along second dimension

                    // reduce dimensions for next level
                    sl1 = sl1 / 2;
                    sl2 = sl2 / 2;

                } // loop over levels of DWT

            } // omp parallel

        }

        void
        unpackdouble	(const T * const x, const int n, const int nc, const int k, T * const y)
        {
            for (int i = 0; i < n; i++)
            {
                y [i] = x [k + nc * i];
            }
        }

        void
        packdouble		(const T * const x, const int n, const int nc, const int k, T * const y)
        {
            for (int i = 0; i < n; i++)
                y [k + nc * i] = x [i];
        }


        void
        copydouble		(const T * const src, T * const dest, const int n)
        {
        	memcpy (dest, src, n * sizeof (T));
        }


        void
        adddouble       (const T * const x, const T * const y, const int n, T * const z)
        {
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
        void
        downhi			(const T * const signal, const int side_length, T * const dwt_high)
        {

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
            for (int i = mlo; i < n2; i++)
            {

                // result of convolution
                s = 0.;

                // perform convolution
                for (int h = 0; h < _fl; h++)
                    s += _hpf_d [h] * signal [2 * i + 1 - h];

                // assign result of convolution
                dwt_high [i] = s;

            } // loop over pixels of dwt_high

            // case: filter_length > side_length => only edge values
            if (mlo > n2)
                mlo = n2;


            /* fix up edge values */

            // loop over edge values
            for (int i = 0; i < mlo; i++)
            {

                // result of convolution
                s = 0.;

                // start signal index for convolution
                j = 2 * i + 1;

                // loop over filter elements
                for (int h = 0; h < _fl; h++)
                {

                    // adjust index if it exceeds side_length
                    if (j < 0)
                        j += side_length;

                    s += _hpf_d [h] * signal [j];

                    // update index
                    --j;

                }

                // assign result of convolution
                dwt_high [i] = s;

            } // loop over edge values

        }


        /**
         * @brief               Perform lowpass part of one level forward 1D DWT on signal with given side_length.
         *
         * @param  signal       Signal to be transformed.
         * @param  side_length  Side length of current DWT level.
         * @param  dwt_low      Resulting DWT.
         */
        void
        downlo  		        (const T * const signal, const int side_length, T * const dwt_low)
        {

            int j;
            T s;

            /*lowpass version */

            // half of side_length (length of dwt_low)
            const int n2 = side_length / 2;

            // half of filter_length
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
            for (int i= 0; i <= mhi; i++)
            {

                // result of convolution
                s = 0.;
                // apply low pass filter (convolution)
                for (int h = 0; h < _fl; h++)
                    s += _lpf_d [h] * signal [2 * i + h];
                // assign result of convolution
                dwt_low [i] = s;

            } // loop over pixels of dwt_low


            /* fix up edge values */

            // loop over edge values (periodic boundary)
            for (int i = mhi + 1; i < n2; i++)
            {

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


        void
        idpwt2		(const Matrix <T> & wc, const int ell, const int J,
              		 container<T>& temp_mem, Matrix <T> & img)
        {

            T * wcplo, * wcphi, * templo, * temphi, * temptop, * tmp;

            // assign dwt to result image
            img = wc;

            // calculate start level for backwards DWT
            int sl1 = _sl1_scale,
                sl2 = _sl2_scale,
                sl3 = _sl3_scale;

# pragma omp parallel default (shared) firstprivate (sl1, sl2, sl3) \
                     private (wcplo, wcphi, temphi, templo, temptop, tmp) num_threads (_num_threads)
            {

                size_t stride;
                const int t_num = omp_get_thread_num ();

                // loop over levels of backwards DWT
                for (int j = ell; j < J; j++)
                {

                    // update stride
                    stride = 6 * sl2 * t_num;
                    tmp = & temp_mem [stride];
                    templo  = & temp_mem [2 * sl2 + stride];
                    temphi  = & temp_mem [3 * sl2 + stride];
                    temptop = & temp_mem [4 * sl2 + stride];

# pragma omp for schedule (OMP_SCHEDULE)
                    // loop over lines along second dimension ('rows') of result image
                    for (int c1 = 0; c1 < 2 * sl1; c1++)
                    {

                        // copy lowpass part of current line to temporary memory
                        unpackdouble (& img [0], sl2, _sl1, c1, templo);

                        // copy highpass part of current line to temporary memory
                        unpackdouble (& img [sl2 * _sl1], sl2, _sl1, c1, temphi);

                        // perform lowpass reconstruction
                        uplo (templo, sl2, tmp);
                        // perform highpass reconstruction
                        uphi (temphi, sl2, temptop);

                        // fusion of reconstruction parts
                        adddouble (tmp, temptop, sl2 * 2, tmp);

                        // write back reconstructed line
                        packdouble (tmp, sl2 * 2, _sl1, c1, & img [0]);

                    } // loop over lines along second dimension of result image

                    // update stride
                    stride = 5 * sl1 * t_num;
                    tmp = & temp_mem [stride];
                    templo = & temp_mem [    sl1 + stride];
                    temphi = & temp_mem [3 * sl1 + stride];

# pragma omp for schedule (OMP_SCHEDULE)
                    // loop  over lines along first dimension ('columns') of result image
                    for (int c2 = 0; c2 < 2 * sl2; c2++)
                    {

                        // assign address of current line's lowpass part
                        wcplo = & img [c2 * _sl1];
                        // assign address of current line's highpass part
                        wcphi = & img [c2 * _sl1 + sl1];

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


        void
        uplo		(const T * const wc, const int side_length, T * const signal)
        {

            int j;
            T s, s_odd;

            /*lowpass version */

            /* away from edges */

            // upper bound for even filter indices
            const int meven = (_fl + 1) / 2;
            // upper bound for odd filter indices
            const int modd = _fl / 2;

            // loop over regular signal indices
            for (int i = meven; i < side_length; i++)
            {

                // init convolution results
                s = 0.;
                s_odd = 0.;

                // perform convolution for even and odd filter indices
                for (int h = 0; h < modd; h++)
                {

                    // even filter index
                    s += _lpf_r [2 * h] * wc [i - h];
                    // odd filter index
                    s_odd += _lpf_r [2 * h + 1] * wc [i - h];

                }
                // case of odd filter_length (-> more even indices: start with index 0)
                if (meven > modd)
                    s += _lpf_r [2 * meven] * wc [i - meven];

                // assign convolution results
                signal [2 * i] = s;
                signal [2 * i + 1] = s_odd;

            } // loop over regular signal indices


            /* fix up edge values */

            // upper bound for filter indices
            int mmax = meven;
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
                for (int h = 0; h < modd; h++)
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
                if (meven > modd)
                    s += _lpf_r [2 * meven] * wc [j];

                // assign convolution results
                signal [2 * i] = s;
                signal [2 * i + 1] = s_odd;

            } // loop over edge values

        }


        void
        uphi		(const T * const wc, const int side_length, T * const signal)
        {

            int j;
            T s, s_odd;

            /*hipass version */

            // upper bound for even filter indices
            const int meven = (_fl + 1) / 2;
            // upper bound for odd filter indices
            const int modd = _fl / 2;

            /* away from edges */

            // loop over regular signal indices
            for (int i = 0; i < side_length - meven; i++)
            {

                // init convolution results
                s = 0.;
                s_odd = 0.;

                // perform convolution for even and odd filter indices
                for (int h = 0; h < modd; h++)
                {

                    // even filter index
                    s += _hpf_r [2 * h] * wc [i + h];
                    // odd filter index
                    s_odd += _hpf_r [2 * h + 1] * wc [i + h];

                } // perform convolution

                // case of odd filter_length
                if (meven > modd)
                    s += _hpf_r [2 * meven] * wc [i + meven];

                // assign convolution results
                signal [2 * i + 1] = s;
                signal [2 * i] = s_odd;

            } // loop over regular signal indices


            /* fix up edge values */

            // lower bound for indices of edge values
            int mmin = side_length - meven;
            // possible correction if mmin less than zero
            if (mmin < 0)
                mmin = 0;

            // loop over edge values
            for (int i = mmin; i < side_length; i++)
            {

                // init convolution results
                s = 0.;
                s_odd = 0.;
                // start index of wavelet coefficients
                j = i;

                // perform convolution for even and odd indices
                for (int h = 0; h < meven; h++)
                {

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
