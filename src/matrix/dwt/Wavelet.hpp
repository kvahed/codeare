# ifndef __WAVELET_HPP__


# define __WAVELET_HPP__


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
template <> template <class T, int wl_mem>
class WaveletTraits <T, WL_HAAR, wl_mem>
{

    public:

        /**
         * @brief               Getter for low pass reconstruction filter of haar wavelet.
         *
         * @return              Low pass reconstruction filter.
         */
        static inline
        void
        LowPassFilterDecom      (Matrix <T> & lpf_d)
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
        static inline
        void
        HighPassFilterDecom     (Matrix <T> & hpf_d)
        {
            Matrix <T> lpf_d (2);
            LowPassFilterDecom (lpf_d);
            mirrorfilt (lpf_d, hpf_d);
        }

        /**
         * @brief               Getter for low pass decomposition filter of haar wavelet.
         *
         * @return              Low pass decomposition filter.
         */
        static inline
        void
        LowPassFilterRecon      (Matrix <T> & lpf_r)
        {
            LowPassFilterDecom (lpf_r);
        }

        /**
         * @brief               Getter for high pass decomposition filter of haar wavelet.
         *
         * @return              High pass decomposition filter.
         */
        static inline
        void
        HighPassFilterRecon     (Matrix <T> & hpf_r)
        {
            HighPassFilterDecom (hpf_r);
        }

};


/**
 * @brief           Implementation of Daubechies wavelet, member: 4.
 */
template <>
template <class T>
class WaveletTraits <T, WL_DAUBECHIES, 4>
{

    public:

        /**
         * @brief               Getter for low pass reconstruction filter of Daubechies wavelet, member: 4.
         *
         * @return              Low pass reconstruction filter.
         */
        static inline
        void
        LowPassFilterDecom      (Matrix <T> & lpf_d)
        {
            lpf_d [0] =  0.48296291314453414337487159986;
            lpf_d [1] =  0.83651630373780790557529378092;
            lpf_d [2] =  0.22414386804201338102597276224;
            lpf_d [3] = -0.12940952255126038117444941881;
        }

        /**
         * @brief               Getter for high pass reconstruction filter of Daubechies wavelet, member: 4.
         *
         * @return              High pass reconstruction filter.
         */
        static inline
        void
        HighPassFilterDecom     (Matrix <T> & hpf_d)
        {
            Matrix <T> lpf_d (4);
            LowPassFilterDecom (lpf_d);
            mirrorfilt (lpf_d, hpf_d);
        }

        /**
         * @brief               Getter for low pass decomposition filter of Daubechies wavelet, member: 4.
         *
         * @return              Low pass decomposition filter.
         */
        static inline
        void
        LowPassFilterRecon      (Matrix <T> & lpf_r)
        {
            LowPassFilterDecom (lpf_r);
        }

        /**
         * @brief               Getter for high pass decomposition filter of Daubechies wavelet, member: 4.
         *
         * @return              High pass decomposition filter.
         */
        static inline
        void
        HighPassFilterRecon     (Matrix <T> & hpf_r)
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

        /**
         * @brief               Getter for low pass reconstruction filter of Daubechies wavelet, member: 8.
         *
         * @return              Low pass reconstruction filter.
         */
        static inline
        void
        LowPassFilterDecom      (Matrix <T> & lpf_d)
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
        static inline
        void
        HighPassFilterDecom     (Matrix <T> & hpf_d)
        {
            Matrix <T> lpf_d (8);
            LowPassFilterDecom (lpf_d);
            mirrorfilt (lpf_d, hpf_d);
        }

        /**
         * @brief               Getter for low pass decomposition filter of Daubechies wavelet, member: 8.
         *
         * @return              Low pass decomposition filter.
         */
        static inline
        void
        LowPassFilterRecon      (Matrix <T> & lpf_r)
        {
            LowPassFilterDecom (lpf_r);
        }

        /**
         * @brief               Getter for high pass decomposition filter of Daubechies wavelet, member: 8.
         *
         * @return              High pass decomposition filter.
         */
        static inline void
        HighPassFilterRecon     (Matrix <T> & hpf_r) {
            HighPassFilterDecom (hpf_r);
        }

};


/**
 * @brief               Construct mirror filter of given low pass filter.
 */
template <class T>
static inline
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


/**
 * @brief				Setup given filter matrices for specified wavlet.
 */
template <class T>
static
void
setupWlFilters			(const wlfamily wl_fam, const int wl_mem, Matrix <T> & lpf_d, Matrix <T> & lpf_r,
							Matrix <T> & hpf_d, Matrix <T> & hpf_r)
{
	switch (wl_fam)
	{
	case WL_DAUBECHIES:
		switch (wl_mem)
		{
		case 8:
            WaveletTraits <T, WL_DAUBECHIES, 8> :: LowPassFilterDecom (lpf_d);
            WaveletTraits <T, WL_DAUBECHIES, 8> :: LowPassFilterRecon (lpf_r);
            WaveletTraits <T, WL_DAUBECHIES, 8> :: HighPassFilterDecom (hpf_d);
            WaveletTraits <T, WL_DAUBECHIES, 8> :: HighPassFilterRecon (hpf_r);
            break;
		case 4:
			WaveletTraits <T, WL_DAUBECHIES, 4> :: LowPassFilterDecom (lpf_d);
			WaveletTraits <T, WL_DAUBECHIES, 4> :: LowPassFilterRecon (lpf_r);
			WaveletTraits <T, WL_DAUBECHIES, 4> :: HighPassFilterDecom (hpf_d);
			WaveletTraits <T, WL_DAUBECHIES, 4> :: HighPassFilterRecon (hpf_r);
			break;
		}
		break;
	case WL_HAAR:
        WaveletTraits <T, WL_HAAR, 2> :: LowPassFilterDecom (lpf_d);
        WaveletTraits <T, WL_HAAR, 2> :: LowPassFilterRecon (lpf_r);
        WaveletTraits <T, WL_HAAR, 2> :: HighPassFilterDecom (hpf_d);
        WaveletTraits <T, WL_HAAR, 2> :: HighPassFilterRecon (hpf_r);
        break;
	}
}


# endif // __WAVELET_HPP__
