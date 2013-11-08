# ifndef __OCL_DWT_HPP__

# define __OCL_DWT_HPP__


/**
 * includes
 */
# include "oclmatrix/oclMatrix.hpp"
# include "oclTraits.hpp"
# include "oclSettings.hpp"

/**
 * @brief 3D Discrete wavelet transform for oclMatrix
 */
template <class T>
class oclDWT
{


    public:

        /**
         * @brief           Construct 3D Wavelet transform
         *                  with wavelet class and side length
         *
         * @param  wf       Wavelet family (default none, i.e. ID)
         * @param  wm       Familty member (default 4)
         */
        oclDWT (const wlfamily & wf = ID,
                const   size_t & wm = 4)
        {

            print_optional ("oclDWT :: oclDWT (...)", v_level);

            this->m_wf = wf;

        }


        virtual
        ~oclDWT ()
        {

            print_optional ("oclDWT :: ~oclDWT ()", v_level);

        }


        /**
         * @brief       Forward transform
         *
         * @param  m    To transform
         *
         * @return      Transform
         */
        oclMatrix <T>
        Trafo           (const oclMatrix <T> & m)
        const
        {

            // create result matrix
            oclMatrix <T> res (m.Height (), m.Width ());

            // calculate dwt
            oclOperations <T> :: ocl_function_dwt_forward (m.mp_oclData, res.mp_oclData, m.Height (), m.Width ());

            // return result
            return res;

        }


        /**
         * @brief       Adjoint transform
         *
         * @param  m    To transform
         *
         * @return      Transform
         */
        oclMatrix <T>
        Adjoint         (const oclMatrix <T> & m)
        const
        {

            /* oclOperation */

        }


        /**
         * @brief       Forward transform
         *
         * @param  m    To transform
         *
         * @return      Transform
         */
        inline
        oclMatrix <T>
        operator*       (const oclMatrix <T> & m)
        const
        {

            return Trafo (m);

        }


        /**
         * @brief       Adjoint transform
         *
         * @param  m    To transform
         *
         * @return      Transform
         */
        inline
        oclMatrix <T>
        operator->* (const oclMatrix <T> & m)
        const
        {

            return Adjoint (m);

        }



    private:

        wlfamily m_wf; /**< @brief wavelet family */

        size_t m_sz; /**< @brief data size */
        size_t m_sl; /**< @brief side length */

        /* private member for verbosity level of class */
        static const VerbosityLevel v_level;

};



/*************************************
 ** initialize static class members **
 *************************************/
template <class T>
const VerbosityLevel oclDWT <T> :: v_level = global_verbosity [OCL_DWT];


# endif __OCL_DWT_HPP__
