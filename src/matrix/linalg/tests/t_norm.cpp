#include "Matrix.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Lapack.hpp"
#include "Print.hpp"

#define VERBOSE

template<class T> bool norm_check () {

    Matrix<T> A = rand<T> (4,3);
	double nrm1 = norm(A), nrm2 = 0.;

	for (size_t i = 0; i < A.Size(); ++i)
		nrm2 += TypeTraits<T>::Real(A[i]*TypeTraits<T>::Conj(A[i]));
	nrm2 = sqrt(nrm2);
    
#ifdef VERBOSE
    std::cout << "A=\n" << A;
    std::cout << "norm(A)=" << nrm1 << " " << nrm2 << std::endl;
#endif

	return ((nrm1-nrm2)/nrm1 < 1.e-6);

}

int main (int args, char** argv) {

	return  (norm_check<cxfl>()  && norm_check<cxdb>() && 
			 norm_check<float>() && norm_check<double>()) ? 0 : 1;
		
}
