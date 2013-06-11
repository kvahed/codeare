/*
 *  codeare Copyright (C) 2013 Daniel Joergens
 *
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


# ifndef __DWT2_HPP__

# define __DWT2_HPP__


/**
 * @brief  Supported wavelet families
 */


# include "Matrix.hpp"
# include "Wavelet.hpp"



/**
 * @brief 2D Discrete wavelet transform for Matrix template (from GSL)
 */
template <   class         T,
          wlfamily    wl_fam = WL_DAUBECHIES,
               int    wl_mem = 8,
             class wl_traits = WaveletTraits <T, wl_fam, wl_mem> >
class DWT2 {


    public:


        /**
         * @brief Construct 2D Wavelet transform with wavelet class and side length
         *
         * @param  sl      Side length
         */
        DWT2 (const size_t sl, const int wl_scale = 0, const int dim = 2)
            : m_dim (dim),
              m_lpf_d (wl_mem),
              m_lpf_r (wl_mem),
              m_hpf_d (wl_mem),
              m_hpf_r (wl_mem),
              m_max_size (sl * sl),
              temp ((T *) malloc (4 * m_max_size * sizeof (T))),
              m_wl_scale (wl_scale)
        {

            wl_traits::LowPassFilterDecom (m_lpf_d);
            wl_traits::LowPassFilterRecon (m_lpf_r);
            wl_traits::HighPassFilterDecom (m_hpf_d);
            wl_traits::HighPassFilterRecon (m_hpf_r);

            omp_set_num_threads (2);

        }


        virtual
        ~DWT2 ()
        {
            /* -- */

            free (temp);

        }


        /**
         * @brief    Forward transform
         *
         * @param  m To transform
         * @return   Transform
         */
        inline
        Matrix<T>
        Trafo        (const Matrix<T>& m)
        {
            return TransformForward (m);
        }


        /**
         * @brief    Adjoint transform
         *
         * @param  m To transform
         * @return   Transform
         */
        inline
        Matrix<T>
        Adjoint      (const Matrix<T>& m)
        {
            return TransformBackwards (m);
        }


        /**
         * @brief    Forward transform
         *
         * @param  m To transform
         * @return   Transform
         */
        inline
        Matrix<T>
        operator*    (const Matrix<T>& m)
        {
            return Trafo(m);
        }


        /**
         * @brief    Adjoint transform
         *
         * @param  m To transform
         * @return   Transform
         */
        inline
        Matrix<T>
        operator->* (const Matrix<T>& m)
        {
            return Adjoint(m);
        }


    private:


        /**
         * type definitions
         */
        typedef typename elem_type_traits <T> :: value_type value_type;


        /**
         * variable definitions
         */

        // dimension of DWT2
        const int m_dim;

        // low pass filters
        Matrix <value_type> m_lpf_d;
        Matrix <value_type> m_lpf_r;

        // high pass filters
        Matrix <value_type> m_hpf_d;
        Matrix <value_type> m_hpf_r;

        // maximum size of matrix
        const size_t m_max_size;

        // temporary memory
        T * temp;

        // wavelet scale (max. decomposition level)
        const int m_wl_scale;


        /**
         * function definitions
         */

        /**
         * @brief   Transform
         *
         * @param   m   To transform
         * @param   bw  Backward: true, Forward: false
         * @return      Transform
         */
        inline
        Matrix<T>
        TransformForward    	(const Matrix<T>& m)
        {

            assert (m.Size () <= m_max_size);

            Matrix<T> res (m.DimVector());

//            if (m_dim == 2)
//            {

                // create vars from mex function
                int J = 0, nn;
                for (nn = 1; nn < m.Height (); nn *= 2 )
                    J ++;
                if (nn  !=  m.Height ()){
                    std::cout << "FWT2 requires dyadic length sides" << std::endl;
                    assert (false);
                }

                // call dpwt2
                res = dpwt2 (m, m_wl_scale, J, temp);

//            }

            return res;

        }


        /**
         * backwards transform
         */
        inline
        Matrix <T>
        TransformBackwards		(const Matrix <T> & m)
        {

            assert (m.Size () <= m_max_size);

            Matrix <T> res (m.DimVector());

//            if (m_dim == 2)
//            {

                // create vars from mex function
                int J = 0, nn;
                for (nn = 1; nn < m.Height (); nn *= 2 )
                    J ++;
                if (nn  !=  m.Height ()){
                    std::cout << "IWT2 requires dyadic length sides" << std::endl;
                    assert (false);
                }

                // call dpwt2
                res = idpwt2 (m, m_wl_scale, J, temp);

//            }

            return res;

        }


        /**
         *  COPIED FROM WAVELAB IMPLEMENTATION
         */

        Matrix <T>
        dpwt2		(const Matrix <T> & sig, const int ell, const int J, T * temp)
        {

            Matrix <T> res (sig.DimVector ());

            T *wcplo,*wcphi,*templo,*temphi;
            //		copydouble(&sig[0],&res[0],sig.Height()*sig.Width());

            // assign signal to result matrix
            res = sig;

            templo = &temp[sig.Height ()];
            temphi = &temp[2*sig.Height()];

            int num_rows = sig.Height();
            int num_cols = sig.Width();
            int side_length = sig.Height();

            // loop over levels of DWT
            for (int j = (J-1); j >= ell; --j)
            {

                // loop over columns of image
                for (int col=0; col < side_length; col++)
                {

                    // access to lowpass part of DWT
                    wcplo = &res[col*num_rows];
                    // access to highpass part of DWT
                    wcphi = &res[col*num_rows + side_length/2];

                    // copy part of image to temp memory
                    copydouble (wcplo, temp, side_length);

                    // apply low pass filter on column and write to result matrix
                    downlo (temp, side_length, wcplo);
                    // apply high pass filter on column and wirte to result matrix
                    downhi (temp, side_length, wcphi);

                } // loop over columns

                // loop over rows of image
                for (int row=0; row < side_length; row++)
                {

                    // copy row-th row of imag to temp
                    unpackdouble (res.Memory(0), side_length, num_cols, row, temp);

                    // apply low pass filter on row and write to temp mem
                    downlo (temp, side_length, templo);
                    // apply high pass filter on row and write to temp mem
                    downhi (temp, side_length, temphi);

                    // write temp lowpass result to result matrix
                    packdouble (templo, side_length/2, num_cols, row, &res[0]);
                    // write temp highpass result to result matrix
                    packdouble (temphi, side_length/2, num_cols, row, &res[side_length/2*num_rows]);

                } // loop over rows of image

                // reduce dimension for next level
                side_length = side_length/2;

            } // loop over levels of DWT

            return res;

        }

        void
        unpackdouble	(const T * const x, const int n, const int nc, const int k, T * const y)
        {

            int i;

//# pragma omp parallel for num_threads (2)
            for( i=0; i < n; i++)
                y[i] = x[k+nc*i];
        }

        void
        packdouble		(const T * const x, const int n, const int nc, const int k, T * const y)
        {
            int i;

//# pragma omp parallel for num_threads (2)
            for( i=0; i < n; i++)
                y[k+nc*i] = x[i];
        }


        void
        copydouble		(const T * const src, T * const dest, const int n)
        {

            // TODO: use operator= instead ...

//# pragma omp parallel for num_threads (2)
            for (int i = 0; i < n; ++i)
            {
                dest [i] = src [i];
            }

        }


        void
        adddouble       (const T * const x, const T * const y, const int n, T * const z)
        {

            // TODO: use operator+ instead ...
//# pragma omp parallel for num_threads (2)
            for (int i = 0; i < n; ++i)
            {
                z[i] = x[i] + y[i];
            }

        }


        /**
         * @brief               Perform highpass part of one level forward 1D DWT on signal with given side_length.
         *
         * @param  signal       Signal to be transformed.
         * @param  side_length  Side length of current DWT level.
         * @param  dwt_hgih     Resulting DWT.
         */
        void
        downhi			(const T * const signal, const int side_length, T * const dwt_high)
        {
            int n2, mlo, j;
            T s;

            int filter_length = m_hpf_d.Dim (0);

            /* highpass version */

            // half of side_length
            n2 = side_length/2;

            // half of filter length
            mlo = filter_length/2 - 1;

            // adjust lower bound if to low
            if (2*mlo+1 - (filter_length-1) < 0)
                mlo++;

            // loop over pixels of dwt_high
//# pragma omp parallel for private (s) num_threads (2)
            for (int i = mlo; i < n2; i++)
            {

                // result of convolution
                s = 0.;

                // perform convolution
                for (int h = 0; h < filter_length; h++)
                    s += m_hpf_d [h]* signal [2*i+1-h];

                // assign result of convolution
                dwt_high[i] = s;

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
                j = 2*i+1;

                // loop over filter elements
                for (int h = 0; h < filter_length; h++)
                {

                    // adjust index if it exceeds side_length
                    if (j < 0)
                        j += side_length;

                    s += m_hpf_d [h] * signal [j];

                    // update index
                    --j;

                }

                // assign result of convolution
                dwt_high[i] = s;

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

            int n2, mlo, mhi, j;
            T s;

            int filter_length = m_lpf_d.Dim (0);

            /*lowpass version */

            // half of side_length (length of dwt_low)
            n2 = side_length/2;

            // half of filter_length
            mlo = filter_length /2;

            // upper bound for "normal" convolution
            mhi = n2 - mlo;

            // upper bound to high
            if (2*mhi + (filter_length-1) >= side_length)
                --mhi;
            // upper bound to low
            if (mhi < 0)
                mhi = -1;

            // loop over pixels of dwt_low
//# pragma omp parallel for private (s) num_threads (2)
            for (int i= 0; i<=mhi; i++)
            {

                // result of convolution
                s = 0.;
                // apply low pass filter (convolution)
                for (int h = 0; h < filter_length; h++)
                    s += m_lpf_d [h] * signal [2*i+h];
                // assign result of convolution
                dwt_low [i] = s;

            } // loop over pixels of dwt_low


            /* fix up edge values */

            // loop over edge values (periodic boundary)
            for (int i = mhi+1; i < n2; i++)
            {

                // result of convolution
                s = 0.;

                // start signal index for convolution
                j = 2*i;

                // loop over filter elements
                for (int h = 0; h < filter_length; h++){

                    // adjust index if it exceeds current side_length
                    if (j >= side_length)
                        j -= side_length;
                    s += m_lpf_d [h] * signal [j];

                    // update index
                    j++;

                } // perform convolution

                // assign result of convolution
                dwt_low [i] = s;

            }

        }


        Matrix <T>
        idpwt2		(const Matrix <T> & wc, const int ell, const int J, T * temp)
        {

            Matrix <T> img (wc.DimVector());

            const int nr = wc.Height();
            const int nc = wc.Width();

            T *wcplo,*wcphi,*templo,*temphi,*temptop;
            int nj;
            //		copydouble(wc,img,nr*nc);

            // assign dwt to result image
            img = wc;

            templo = &temp[nr];
            temphi = &temp[2*nr];
            temptop = &temp[3*nr];

            // calculate start level for backwards DWT
            nj = 1;
            for (int k = 0; k < ell; k++)
                nj *=2;

            // loop over levels of backwards DWT
            for (int j = ell; j < J; j++)
            {

                // loop over rows of result image
                for (int k = 0; k < 2 * nj; k++)
                {

                    // copy lowpass part of current row to temporary memory
                    unpackdouble(&img[0],nj,nc,k,templo);
                    // copy highpass part of current row to temporary memory
                    unpackdouble(&img[nj*nr],nj,nc,k,temphi);

                    // perform lowpass reconstruction
                    uplo(templo, nj,temp);
                    // perform highpass reconstruction
                    uphi(temphi, nj, temptop);

                    // fusion of reconstruction parts
                    adddouble(temp,temptop,nj*2,temp);

                    // write back reconstructed row
                    packdouble(temp,nj*2,nc,k,&img[0]);

                } // loop over rows of result image

                // loop  over cols of result image
                for (int k = 0; k < 2 * nj; k++)
                {

                    // assign address of current column's lowpass part
                    wcplo = &img[k*nr];
                    // assign address of current column's highpass part
                    wcphi = &img[k*nr + nj];

                    // copy lowpass part to temporary memory
                    copydouble(wcplo,temp,nj);

                    // perform lowpass reconstruction
                    uplo(wcplo, nj, templo);
                    // perform highpass reconstruction
                    uphi(wcphi, nj, temphi);

                    // combine reconstructed parts and write back to current column
                    adddouble(templo,temphi,nj*2,wcplo);

                } // loop over cols of result image

                // update current row / column size
                nj *= 2;

            } // loop over levels of backwards DWT

            return img;

        }


        void
        uplo		(const T * const wc, const int side_length, T * const signal)
        {
            int j, meven, modd, mmax;
            T s, s_odd;

            const int filter_length = m_lpf_r.Dim (0);

            /*lowpass version */

            /* away from edges */

            // upper bound for even filter indices
            meven = (filter_length+1)/2;
            // upper bound for odd filter indices
            modd = filter_length/2;

            // loop over regular signal indices
//# pragma omp parallel for private (s, s_odd) num_threads (2)
            for (int i = meven; i < side_length; i++)
            {

                // init convolution results
                s = 0.;
                s_odd = 0.;

                // perform convolution for even and odd filter indices
                for (int h = 0; h < modd; h++)
                {

                    // even filter index
                    s += m_lpf_r [2*h] * wc [i-h];
                    // odd filter index
                    s_odd += m_lpf_r [2*h+1] * wc [i-h];

                }
                // case of odd filter_length (-> more even indices: start with index 0)
                if (meven > modd)
                    s += m_lpf_r [2*meven] * wc [i-meven];

                // assign convolution results
                signal [2*i] = s;
                signal [2*i+1] = s_odd;

            } // loop over regular signal indices


            /* fix up edge values */

            // upper bound for filter indices
            mmax = meven;
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
                    s += m_lpf_r [2*h] * wc [j];
                    // odd part of convolution
                    s_odd += m_lpf_r [2*h+1] * wc [j];

                    // update index
                    --j;

                } // perform convolution

                // case of odd filter_length
                if (meven > modd)
                    s += m_lpf_r [2*meven] * wc [j];

                // assign convolution results
                signal [2*i] = s;
                signal [2*i+1] = s_odd;

            } // loop over edge values

        }


        void
        uphi		(const T * const wc, const int side_length, T * const signal)
        {
            int  meven, modd, j, mmin;
            T s, s_odd;

            const int filter_length = m_hpf_r.Dim (0);

            /*hipass version */

            // upper bound for even filter indices
            meven = (filter_length+1)/2;
            // upper bound for odd filter indices
            modd = filter_length/2;

            /* away from edges */

            // loop over regular signal indices
//# pragma omp parallel for private (s, s_odd) num_threads (2)
            for (int i = 0; i < side_length - meven; i++)
            {

                // init convolution results
                s = 0.;
                s_odd = 0.;

                // perform convolution for even and odd filter indices
                for (int h = 0; h < modd; h++)
                {

                    // even filter index
                    s += m_hpf_r [2*h] * wc [i+h];
                    // odd filter index
                    s_odd += m_hpf_r [2*h+1] * wc [i+h];

                } // perform convolution

                // case of odd filter_length
                if (meven > modd)
                    s += m_hpf_r [2*meven] * wc [i+meven];

                // assign convolution results
                signal [2*i+1] = s;
                signal [2*i] = s_odd;

            } // loop over regular signal indices


            /* fix up edge values */

            // lower bound for indices of edge values
            mmin = side_length - meven;
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
                    s += m_hpf_r [2*h] * wc [j];
                    // odd filter index
                    s_odd += m_hpf_r [2*h+1] * wc [j];

                    // update current wc index
                    j++;

                } // perform convolution

                // assign convolution results
                signal [2*i+1] = s;
                signal [2*i] = s_odd;

            } // loop over edge values

        }



};


# endif // __DWT2_HPP__
