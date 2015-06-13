/*
 * CGLS.cpp
 *
 *  Created on: Apr 5, 2015
 *      Author: kvahed
 */
#include "CGLS.hpp"

namespace codeare {
    namespace optimisation {

template class CGLS<float>;
template class CGLS<double>;
template class CGLS<std::complex<float> >;
template class CGLS<std::complex<double> >;

    }}




