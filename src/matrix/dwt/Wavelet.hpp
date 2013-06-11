# ifndef __WAVELET_HPP__

/************
 ** makros **
 ************/
# define __WAVELET_HPP__


/**************
 ** includes **
 **************/
# include "ElemTypeTraits.hpp"
# include "DWT.hpp"


/**
 * @brief           Templated base class for supported wavelets.
 */
template <class T, wlfamily wl_fam, int wl_mem>
class WaveletTraits
{ /* -- */ };


/**
 * @brief           Implementation of haar wavelet.
 */
template <>
template <class T, int wl_mem>
class WaveletTraits <T, WL_HAAR, wl_mem>
{

    public:

        // extract "real" value type
        typedef typename elem_type_traits <T> :: value_type value_type;

        /**
         * @brief               Getter for low pass reconstruction filter of haar wavelet.
         *
         * @return              Low pass reconstruction filter.
         */
        static
        inline
        void
        LowPassFilterDecom      (Matrix <value_type> & lpf_d)
        {
            float norm_factor = 1 / sqrt (2);

            lpf_d [0] = 1; lpf_d [1] = 1;

            lpf_d /= norm_factor;

        }

        /**
         * @brief               Getter for high pass reconstruction filter of haar wavelet.
         *
         * @return              High pass reconstruction filter.
         */
        static
        inline
        void
        HighPassFilterDecom     (Matrix <value_type> & hpf_d)
        {
            Matrix <value_type> lpf_d;
            LowPassFilterDecom (lpf_d);
            mirrorfilt (lpf_d, hpf_d);
        }

        /**
         * @brief               Getter for low pass decomposition filter of haar wavelet.
         *
         * @return              Low pass decomposition filter.
         */
        static
        inline
        void
        LowPassFilterRecon      (Matrix <value_type> & lpf_r)
        {
            LowPassFilterDecom (lpf_r);
        }

        /**
         * @brief               Getter for high pass decomposition filter of haar wavelet.
         *
         * @return              High pass decomposition filter.
         */
        static
        inline
        void
        HighPassFilterRecon     (Matrix <value_type> & hpf_r)
        {
            HighPassFilterDecom (hpf_r);
        }

};


/**
 * @brief           Implementation of Daubechies wavelet, member: 8.
 */
template <>
template <class T>
class WaveletTraits <T, WL_DAUBECHIES, 8>
{

    public:

        // extract "real" value type
        typedef typename elem_type_traits <T> :: value_type value_type;

        /**
         * @brief               Getter for low pass reconstruction filter of Daubechies wavelet, member: 8.
         *
         * @return              Low pass reconstruction filter.
         */
        static
        inline
        void
        LowPassFilterDecom      (Matrix <value_type> & lpf_d)
        {
            lpf_d [0] =  0.23037781330889650086329118304; lpf_d [1] =  0.71484657055291564708992195527;
            lpf_d [2] =  0.63088076792985890788171633830; lpf_d [3] = -0.02798376941685985421141374718;
            lpf_d [4] = -0.18703481171909308407957067279; lpf_d [5] =  0.03084138183556076362721936253;
            lpf_d [6] =  0.03288301166688519973540751355; lpf_d [7] = -0.01059740178506903210488320852;
        }

        /**
         * @brief               Getter for high pass reconstruction filter of Daubechies wavelet, member: 8.
         *
         * @return              High pass reconstruction filter.
         */
        static
        inline
        void
        HighPassFilterDecom     (Matrix <value_type> & hpf_d)
        {
            Matrix <value_type> lpf_d (8);
            LowPassFilterDecom (lpf_d);
            mirrorfilt (lpf_d, hpf_d);
        }

        /**
         * @brief               Getter for low pass decomposition filter of Daubechies wavelet, member: 8.
         *
         * @return              Low pass decomposition filter.
         */
        static
        inline
        void
        LowPassFilterRecon      (Matrix <value_type> & lpf_r)
        {
            LowPassFilterDecom (lpf_r);
        }

        /**
         * @brief               Getter for high pass decomposition filter of Daubechies wavelet, member: 8.
         *
         * @return              High pass decomposition filter.
         */
        static
        inline
        void
        HighPassFilterRecon     (Matrix <value_type> & hpf_r)
        {
            HighPassFilterDecom (hpf_r);
        }

};


/**
 * @brief               Construct mirror filter of given low pass filter.
 */
template <class T>
static
inline
void
mirrorfilt              (const Matrix <T> & lpf, Matrix <T> & hpf)
{
    int isign = 1;
    for (int i = 0; i < lpf.Dim (0); i++)
    {
        hpf [i] = isign * lpf [i];
        isign *= -1;
    }
}


# endif // __WAVELET_HPP__
