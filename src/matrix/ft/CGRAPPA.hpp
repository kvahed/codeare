/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                                  Forschungszentrum Juelich, Germany
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

#ifndef __CGRAPPA_HPP__
#define __CGRAPPA_HPP__

#include "DFT.hpp"
#include "Workspace.hpp"
#include "GRAPPATraits.hpp"
#include "Lapack.hpp"
#include "Print.hpp"
#include "Access.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "IOContext.hpp"


template<class T> inline static bool eq (const Matrix<T>& A, const Matrix<T>& B) NOEXCEPT {
    assert(A.Dim() == B.Dim());
    return A.Container() == B.Container();
}

template<class T> inline static Vector<size_t> _find (const Matrix<T>& M) NOEXCEPT {
    Vector<size_t> ret;
    T zero_t = (T)0;
    for (size_t i = 0; i < M.Size(); ++i)
        if (M[i]!=zero_t)
            ret.push_back(i);
    return ret;
}

template<class T> inline static Matrix<T> col (const Matrix<T>& M) NOEXCEPT {
    return resize(M,prod(size(M)),1);
}

template<class T> inline static T ssum (const Matrix<T>& M) NOEXCEPT {
    return sum(M.Container());
}

/**
 * @brief GRAPPA operator<br/>
 *        Griswold et al. MRM 2002, vol. 47 (6) pp. 1202-1210
 */
template <class T>
class CGRAPPA : public FT<T> {

    typedef typename TypeTraits<T>::RT RT;

public:


    /**
     * @brief          Default constructor
     */
    CGRAPPA() NOEXCEPT :  m_nthreads(1), m_lambda(0), m_nc(1) {}


    /**
     * @brief Construct with parameters
     */
    CGRAPPA (const Params& p) NOEXCEPT {

// Kernel size
        if (p.exists("kernel_size")) {
            try {
                m_kernel = zeros<T>(p.Get<Vector<size_t> >("kernel_size"));
            } catch (const std::exception&) {
                std::cerr << "  WARNING - CGRAPPA: invalid kernel size definition, defaulting to 4x5" << std::endl;
            }
        } else {
            std::cerr << "  WARNING - CGRAPPA: kernel size unspecified, defaulting to 4x5" << std::endl;
            m_kernel = zeros<T>(5,5);
        }
        std::cout << "  kernel size: " << size(m_kernel) << std::endl;

// AC data
        if (p.exists("ac_data")) {
            try {
                m_ac_data = p.Get<Matrix<T> >("ac_data");
            } catch (const std::exception&) {
                std::cerr << "  ERROR - CGRAPPA: auto calibration data is mandatory input.";
                assert (false);
            }
        }
        m_nc = size(m_ac_data,2); // # Coils
        std::cout << "  # coils: " << m_nc << std::endl;

// Tikhonov lambda
        if (p.exists("lambda"))
            m_lambda = fp_cast(p["lambda"]);
        else
            m_lambda = RT(0.);
        std::cout << "  Tikh lambda: " << m_lambda << std::endl;


// Parallelisation
        if (p.exists("nthreads")) {
            m_nthreads = unsigned_cast(p["nthreads"]);
            omp_set_num_threads(m_nthreads);
        } else {
#pragma omp parallel default (shared)
            {
                m_nthreads = omp_get_num_threads();
            }
        }
        std::cout << "  # threads: " << m_nthreads << std::endl;

        CalcCalibMatrix();

    }

    /**
     * @brief    Clean up and destroy
     */
virtual
                                                      ~CGRAPPA () NOEXCEPT {};


    /**
     * @brief    Adjoint transform
     */
    Matrix<T>
    Adjoint (const Matrix<T>& kspace) const NOEXCEPT {

        Matrix<T> res = kspace;
#pragma omp parallel for default (shared)
        for (int coil = 0; coil < (int)m_nc; ++coil)
            Slice(res, coil, ARC(kspace,coil));
        return res;
    }


    /**
     * @brief    Forward transform
     */
    Matrix<T>
    Trafo (const Matrix<T>& image) const NOEXCEPT {
        Matrix<T> res;
        return res;
    }

    /**
     * @brief    Forward transform
     *
     * @param  m To transform
     * @return   Transform
     */
    virtual Matrix<T>
    operator* (const Matrix<T>& m) const NOEXCEPT {
        return Trafo(m);
    }


    /**
     * @brief    Backward transform
     *
     * @param  m To transform
     * @return   Transform
     */
    virtual Matrix<T>
    operator->* (const Matrix<T>& m) const NOEXCEPT {
        return Adjoint (m);
    }

private:

    /**
     * @brief GRAPPA/ARC reconstruction
     *
     * @param  data      Under-sampled measurement reference
     * @param  coil_num  Number of coil to reconstruct for
     * @return           Reconstructed full k-space
     */
    inline Matrix<T> ARC (const Matrix<T>& data, size_t coil_num) const NOEXCEPT {

        size_t max_list_len = 100;
        Vector<size_t> data_size = size(data);
        Vector<size_t> kernel_size = size(m_kernel);
        Matrix<T> under_sampled = zpad (data, data_size[0]+kernel_size[0]-1, data_size[1]+kernel_size[1]-1, m_nc);
        Matrix<T> dummy (kernel_size[0], kernel_size[1], m_nc);
        dummy (kernel_size[0]/2, kernel_size[0]/2, coil_num) = 1.;
        size_t center = _find(dummy)[0];
        Matrix<T> fully_sampled (data_size[0],data_size[1]);
        Matrix<T> kernel, kernels (kernel_size[0]*kernel_size[1]*m_nc,max_list_len);
        Matrix<short> pattern, patterns (kernel_size[0]*kernel_size[1]*m_nc,max_list_len);
        Matrix<T> tmp (kernel_size[0],kernel_size[1],m_nc);

        for (size_t y = 0, list_len = 1; y < data_size[1]; ++y) // Scan k-space for
            for (size_t x = 0, idx = 0; x < data_size[0]; ++x, idx=0) {

                for (size_t ny = 0; ny < kernel_size[1]; ++ny)
                    for (size_t nx = 0; nx < kernel_size[0]; ++nx)
                        for (size_t nc = 0; nc < m_nc; ++nc)
                            tmp (nx,ny,nc) = under_sampled(x+nx,y+ny,nc);
                pattern = col(abs(tmp)>0);

                for (size_t i = 0; i < list_len; ++i)
                    if (pattern.Container()==Column(patterns,i).Container()) {     // Do we know the pattern?
                        idx = i;
                        break;
                    }
                if (idx == 0) {                            // No
                    kernel = Solve (pattern, center);     // Calculate kernel
                    Column (kernels, list_len, kernel);     // Save kernel and pattern
                    Column (patterns, list_len, pattern);
                    list_len++;
                } else
                    kernel = Column (kernels, idx);

                fully_sampled(x,y) = ssum(kernel*col(tmp));

            }

        return fully_sampled;
    }

    // Solve Ax=b
    inline Matrix<T> Solve (Matrix<short> pattern, size_t center) const NOEXCEPT {
        Vector<size_t> kernel_size = size(m_kernel);
        pattern (center) = 0;
        Vector<size_t> pat_ind = _find(pattern);
        Matrix<T> b = m_coil_calib (pat_ind,center);
        Matrix<T> A = m_coil_calib (pat_ind,pat_ind);
        RT lambda = m_lambda*norm(A,'F')/size(A,0);
        Matrix<T> rawkernel = gemm(inv(A + lambda * eye<T>(pat_ind.size())),b);
        Matrix<T> kernel(prod(kernel_size)*m_nc,1);
        for (size_t i = 0; i < pat_ind.size(); ++i)
            kernel[pat_ind[i]] = rawkernel[i];
        return kernel;
    }

    /**
     * @brief Setup calibration matrix
     */
    inline void CalcCalibMatrix () NOEXCEPT {
        Vector<size_t> ac_size = size(m_ac_data);
        Vector<size_t> kernel_size = size(m_kernel);
        Vector<size_t> calib_mat_size (4);
        calib_mat_size[0] = ac_size[0] - kernel_size[0] + 1;
        calib_mat_size[1] = ac_size[1] - kernel_size[1] + 1;
        calib_mat_size[2] = prod (kernel_size);
        calib_mat_size[3] = ac_size[2];
        m_coil_calib = Matrix<T> (calib_mat_size);
        for (size_t j = 0, count = 0; j < kernel_size[1]; ++j)
            for (size_t i = 0; i < kernel_size[0]; ++i, ++count)
                for (size_t m = 0; m < calib_mat_size[1]; ++m)
                    for (size_t l = 0; l < calib_mat_size[0]; ++l)
                        for (size_t n = 0; n < ac_size[2]; ++n)
                            m_coil_calib (l, m, count, n) = m_ac_data (i+l,j+m,n);
        calib_mat_size[0] *= calib_mat_size[1];
        calib_mat_size[1]  = calib_mat_size[2] * calib_mat_size[3];
        calib_mat_size.resize(2);
        m_coil_calib = resize (m_coil_calib, calib_mat_size);
        m_coil_calib = gemm (m_coil_calib, m_coil_calib, 'C');
    }

    Matrix<T>           m_weights; /**< @brief Correction patch     */
    Matrix<T>           m_ac_data; /**< @brief ACS lines            */
    Matrix<T>           m_kernel;  /**< @brief GRAPPA kernel        */
    Matrix<T>           m_coil_calib;

    Matrix<size_t>       m_kdims;   /**< @brief    */
    Matrix<size_t>       m_adims;
    Matrix<size_t>       m_d;       /**< @brief Dimensions           */
    Matrix<size_t>       m_af;      /**< @brief Acceleration factors */
    Matrix<size_t>       m_sdims;   /**< @brief Scan dimensions      */

    RT m_lambda;

    size_t               m_nc;      /**< @brief Number of receive channels */
    size_t               m_nthreads;

};

#endif /* __CGRAPPA_HPP__ */
