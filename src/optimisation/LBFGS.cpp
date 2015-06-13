#include "LBFGS.hpp"

namespace codeare{
    namespace optimisation {

template class LBFGS<float>;
template class LBFGS<double>;
template class LBFGS<std::complex<float> >;
template class LBFGS<std::complex<double> >;

    }}
