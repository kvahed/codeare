/*
 * Triginometry.hpp
 *
 *  Created on: Jan 29, 2013
 *      Author: kvahed
 */
#ifndef TRIGINOMETRY_HPP_
#define TRIGINOMETRY_HPP_

#include "Matrix.hpp"
#include <math.h>

namespace codeare {
    namespace matrix {
        namespace arithmetic {

            template<class T> inline static Matrix<T>
            sin (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::sin(m[i]);
                return ret;
            }

            template<class T> inline static Matrix<T>
            cos (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::cos(m[i]);
                return ret;
            }
            
            template<class T> inline static Matrix<T>
            tan (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::tan(m[i]);
                return ret;
            }

            template<class T> inline static Matrix<T>
            asin (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::asin(m[i]);
                return ret;
            }
            
            template<class T> inline static Matrix<T>
            acos (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::acos(m[i]);
                return ret;
            }
            
            template<class T> inline static Matrix<T>
            atan (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::atan(m[i]);
                return ret;
            }

            template<class T> inline static Matrix<T>
            atan2 (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::atan2(m[i]);
                return ret;
            }

            template<class T> inline static Matrix<T>
            sinh (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::sinh(m[i]);
                return ret;
            }
            
            template<class T> inline static Matrix<T>
            cosh (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::cosh(m[i]);
                return ret;
            }
            
            template<class T> inline static Matrix<T>
            tanh (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::tanh(m[i]);
                return ret;
            }

            template<class T> inline static Matrix<T>
            exp (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::exp(m[i]);
                return ret;
            }
            
            template<class T> inline static Matrix<T>
            log (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::log(m[i]);
                return ret;
            }
            
            template<class T> inline static Matrix<T>
            log10 (const Matrix<T>& m) {
                Matrix<T> ret (m.Dim(), m.Res());
                #pragma parallel for
                for (size_t i = 0; i < numel(m); i++)
                    ret[i] = std::log10(m[i]);
                return ret;
            }
            
        }
    }
}

#endif /* TRIGINOMETRY_HPP_ */
