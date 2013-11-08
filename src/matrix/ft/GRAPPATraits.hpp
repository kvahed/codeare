#ifndef __GRAPPA_TRAITS__
#define __GRAPPA_TRAITS__

template<class T>
struct GRAPPATraits;

extern "C" {

	int s_src_trg_mat (const std::complex<float>*, const int*, const int*, const int&, const int*,
			const std::complex<float>*, const int*, std::complex<float>*, const int*);
	int d_src_trg_mat (const std::complex<double>*, const int*, const int*, const int&, const int*,
			const std::complex<double>*, const int*, std::complex<double>*, const int*);

}


template<>
struct GRAPPATraits<float> {

	typedef float RType;
	typedef std::complex<RType> CType;

	void src_trg_mat (const CType* aln, const int* asz, const int* msz, const int& dim,
		const int* af, const CType* s, const int* ssz, CType* t, const int* tsz) {
		s_src_trg_mat (aln, asz, msz, dim, af, s, ssz, t, tsz);
	}

};

template<>
struct GRAPPATraits<double> {

	typedef double RType;
	typedef std::complex<RType> CType;

	void src_trg_mat (const CType* aln, const int* asz, const int* msz, const int& dim,
		const int* af, const CType* s, const int* ssz, CType* t, const int* tsz) {
		d_src_trg_mat (aln, asz, msz, dim, af, s, ssz, t, tsz);
	}

};

#endif //__GRAPPA_TRAITS__
