# ifndef __WAVELET_HPP__

# define __WAVELET_HPP__


# include "DWT.hpp"


/**
 * @brief               Construct mirror filter of given low pass filter.
 */
template <class T> static inline void
mirrorfilt (T * lpf, T * hpf, const int fl) NOEXCEPT {
    int isign = 1;
    for (int i = 0; i < fl; i++)  {
        hpf [i] = isign * lpf [i];
        isign *= -1;
    }
}


/**
 * @brief           Templated base class for supported wavelets.
 */
template <class T, wlfamily wl_fam, int wl_mem> class WaveletTraits { /* -- */ };


/**
 * @brief           Implementation of haar wavelet.
 */
template <class T, int wl_mem>
class WaveletTraits <T, WL_HAAR, wl_mem> {

    typedef typename TypeTraits <T>::RT RT;

public:

        /**
         * @ brief              Assign memory to wavelet filters.
         */
        static inline void Malloc (RT* & lpf_d, RT* & lpf_r, RT* & hpf_d, RT* & hpf_r) NOEXCEPT {
            RT* tmp = (RT*) malloc (4 * wl_mem * sizeof (typename TypeTraits <T>::RT));
            lpf_d = & tmp [0         ];
            lpf_r = & tmp [    wl_mem];
            hpf_d = & tmp [2 * wl_mem];
            hpf_r = & tmp [3 * wl_mem];
        }

        /**
         * @brief               Getter for low pass reconstruction filter of haar wavelet.
         *
         * @return              Low pass reconstruction filter.
         */
        static inline void DecomFilters (RT* lpf_d, RT* hpf_d)NOEXCEPT {
            RT norm_factor = 1 / sqrt (2.f);
            lpf_d [0] = norm_factor; lpf_d [1] = norm_factor;
            mirrorfilt (lpf_d, hpf_d, 2);
        }

        /**
         * @brief               Getter for low pass decomposition filter of haar wavelet.
         *
         * @return              Low pass decomposition filter.
         */
        static inline void ReconFilters (RT* lpf_r, RT* hpf_r) NOEXCEPT {
            DecomFilters (lpf_r, hpf_r);
        }

};


/**
 * @brief           Implementation of Daubechies wavelet, member: 4.
 */
template <class T> class WaveletTraits <T, WL_DAUBECHIES, 4> {

    private:
        typedef typename TypeTraits <T>::RT RT;

    public:

        /**
         * @ brief              Assign memory to wavelet filters.
         */
        static inline void Malloc  (RT* & lpf_d, RT* & lpf_r, RT* & hpf_d, RT* & hpf_r) NOEXCEPT {
            RT* tmp = (RT*) malloc (4 * 4 * sizeof (typename TypeTraits <T>::RT));
            lpf_d = & tmp [0    ];
            lpf_r = & tmp [    4];
            hpf_d = & tmp [2 * 4];
            hpf_r = & tmp [3 * 4];
        }

        /**
         * @brief               Getter for low pass reconstruction filter of Daubechies wavelet, member: 4.
         *
         * @return              Low pass reconstruction filter.
         */
        static inline void DecomFilters (RT* lpf_d, RT* hpf_d) NOEXCEPT {
            lpf_d[0] = (RT)0.48296291314453414337487159986; lpf_d[1] = (RT) 0.83651630373780790557529378092;
            lpf_d[2] = (RT)0.22414386804201338102597276224; lpf_d[3] = (RT)-0.12940952255126038117444941881;
            mirrorfilt (lpf_d, hpf_d, 4);
        }

        /**
         * @brief               Getter for low pass decomposition filter of Daubechies wavelet, member: 4.
         *
         * @return              Low pass decomposition filter.
         */
        static inline void ReconFilters (RT* lpf_r, RT* hpf_r) NOEXCEPT {
            DecomFilters (lpf_r, hpf_r);
        }

};


/**
 * @brief           Implementation of Daubechies wavelet, member: 8.
 */
template<class T> class WaveletTraits <T, WL_DAUBECHIES, 8> {

    private:
        typedef typename TypeTraits <T>::RT RT;

    public:

        /**
         * @ brief              Assign memory to wavelet filters.
         */
        static inline void Malloc (RT* & lpf_d, RT* & lpf_r, RT* & hpf_d, RT* & hpf_r) NOEXCEPT {
            RT* tmp = (RT*) malloc (4 * 8 * sizeof (typename TypeTraits <T>::RT));
            lpf_d = & tmp [0    ];
            lpf_r = & tmp [    8];
            hpf_d = & tmp [2 * 8];
            hpf_r = & tmp [3 * 8];
        }

        /**
         * @brief               Getter for low pass reconstruction filter of Daubechies wavelet, member: 8.
         *
         * @return              Low pass reconstruction filter.
         */
        static inline void DecomFilters (RT* lpf_d, RT* hpf_d) NOEXCEPT {
            lpf_d[0] = (RT) 0.23037781330889650086329118304; lpf_d[1] = (RT) 0.71484657055291564708992195527;
            lpf_d[2] = (RT) 0.63088076792985890788171633830; lpf_d[3] = (RT)-0.02798376941685985421141374718;
            lpf_d[4] = (RT)-0.18703481171909308407957067279; lpf_d[5] = (RT) 0.03084138183556076362721936253;
            lpf_d[6] = (RT) 0.03288301166688519973540751355; lpf_d[7] = (RT)-0.01059740178506903210488320852;
            mirrorfilt (lpf_d, hpf_d, 8);
        }

        /**
         * @brief               Getter for low pass decomposition filter of Daubechies wavelet, member: 8.
         *
         * @return              Low pass decomposition filter.
         */
        static inline void ReconFilters (RT* lpf_r, RT* hpf_r) NOEXCEPT {
            DecomFilters        (lpf_r, hpf_r);
        }

};


/**
 * @brief				Setup given filter matrices for specified wavlet.
 */
template <class T> static void setupWlFilters (const wlfamily wl_fam, const int wl_mem,
		typename TypeTraits<T>::RT* & lpf_d, typename TypeTraits<T>::RT* & lpf_r,
		typename TypeTraits<T>::RT* & hpf_d, typename TypeTraits<T>::RT* & hpf_r) NOEXCEPT {

	switch (wl_fam) {
    case ID: break;
	case WL_DAUBECHIES:
		switch (wl_mem) {
		case 8:
		    WaveletTraits<T, WL_DAUBECHIES, 8>::Malloc (lpf_d, lpf_r, hpf_d, hpf_r);
            WaveletTraits<T, WL_DAUBECHIES, 8>::DecomFilters (lpf_d, hpf_d);
            WaveletTraits<T, WL_DAUBECHIES, 8>::ReconFilters (lpf_r, hpf_r);
            break;
		case 4:
            WaveletTraits<T, WL_DAUBECHIES, 4>::Malloc (lpf_d, lpf_r, hpf_d, hpf_r);
	        WaveletTraits<T, WL_DAUBECHIES, 4>::DecomFilters (lpf_d, hpf_d);
	        WaveletTraits<T, WL_DAUBECHIES, 4>::ReconFilters (lpf_r, hpf_r);
			break;
		}
		break;
    case WL_DAUBECHIES_CENTERED:
        break;
	case WL_HAAR:
        WaveletTraits<T, WL_HAAR, 2>::Malloc (lpf_d, lpf_r, hpf_d, hpf_r);
        WaveletTraits<T, WL_HAAR, 2>::DecomFilters (lpf_d, hpf_d);
        WaveletTraits<T, WL_HAAR, 2>::ReconFilters (lpf_r, hpf_r);
        break;
    case WL_HAAR_CENTERED:
        break;
    case WL_BSPLINE:
        break;
    case WL_BSPLINE_CENTERED:
        break;
	}
    
}


# endif // __WAVELET_HPP__
