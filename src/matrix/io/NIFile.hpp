/*
 * NIFile.hpp
 *
 *  Created on: May 24, 2013
 *      Author: kvahed
 */

#ifndef NIFILE_HPP_
#define NIFILE_HPP_

#include "Matrix.hpp"
#include "IOFile.hpp"
#include "Algos.hpp"

#include <nifti1_io.h>

template <class T>
struct NITraits;

template<> struct NITraits<float> {
	static const int native = DT_FLOAT;
};
template<> struct NITraits<double> {
	static const int native = DT_DOUBLE;
};
template<> struct NITraits<cxfl> {
	static const int native = DT_COMPLEX;
};
template<> struct NITraits<cxdb> {
	static const int native = DT_COMPLEX128;
};
template<> struct NITraits<short> {
	static const int native = DT_SIGNED_SHORT;
};
template<> struct NITraits<long> {
	static const int native = DT_SIGNED_INT;
};
template<> struct NITraits<size_t> {
	static const int native = (SIZEOF_VOIDP == 8) ? DT_UINT64 : DT_UINT32;
};


using namespace codeare::matrix::io;

class NIFile : public IOFile {



public:


	NIFile (const std::string& fname, const IOMode mode = READ,
			const Params& params = Params(), const bool verbose = false) :
			IOFile (fname, mode, params, verbose) {
		if (mode == READ)
			assert(fexists(fname));
	}


	/**
	 * @brief          Dump matrix to NIFTI file
	 *
	 * @param  M       Matrix
	 * @param  fname   File name
	 * @return         Success
	 */
	template <class T> bool
	Write (const Matrix<T>& M, const std::string& uri = "") {

		Matrix<T> tmp = M;
		tmp = squeeze(tmp);
		size_t l = this->m_fname.length();
		nifti_1_header header;
		header.sizeof_hdr = 348;
		header.dim[0] = ndims(tmp);
		header.pixdim[0] = ndims(tmp);

		if (ndims(tmp) > 7) {
			printf ("Cannot dump more than 8 dimensions to NIFTI file\n.");
			return false;
		}

		for (size_t i = 0; i < 7; ++i) {
			header.dim   [i+1] = tmp.Dim(i);
			header.pixdim[i+1] = tmp.Res(i);
		}

		header.datatype = NITraits<T>::native;

		nifti_image* ni = nifti_convert_nhdr2nim (header, NULL);

		ni->nifti_type = 1;

		// Single nii.gz file
		ni->fname = (char*) calloc(1,l);
		strcpy(ni->fname,this->m_fname.c_str());
		ni->iname = (char*) calloc(1,l);
		strcpy(ni->iname,this->m_fname.c_str());

		ni->data = (void*) malloc (numel(M) * sizeof (T));
		memcpy (ni->data, M.Memory(), numel(M) * sizeof (T));

		nifti_image_write (ni);
		nifti_image_free (ni);

		return true;

	}


	/**
	 * @brief          Read matrix from NIFTI file
	 *
	 * @param  M       Matrix
	 * @param  fname   File name
	 * @return         Success
	 */
	template <class T> Matrix<T>
	Read (const std::string& uri = "") const {

		std::vector<size_t> ndims (INVALID_DIM,1);
		std::vector<float>  nress (INVALID_DIM,1.0);
		size_t i = 0;

		nifti_image* ni = nifti_image_read (this->m_fname.c_str(), 1);
		Matrix<T> M;

		if (ni == NULL)
			return M;

		for (; i < ni->dim[0]; ++i)
			if (ni->dim[i+1] > 1) {
				ndims[i] = ni->dim[i+1];
				nress[i] = ni->pixdim[i+1];
			}

		for (; i < INVALID_DIM; ++i) {
			ndims[i] = 1;
			nress[i] = 1.0;
		}

		M = Matrix<T>(ndims, nress);

		assert (ni->datatype == NITraits<T>::native);
		memcpy (&M[0], ni->data, numel(M)*sizeof(T)); // fast

		nifti_image_free (ni);

		return M;

	}


	/**
	 * @brief Read a particular data set from file
	 *
	 * @return  Success
	 */
	template<class T> Matrix<T>
	Read (const TiXmlElement* txe) const {
		return this->Read<T>("");
	}


	/**
	 * @brief  Write data to file
	 *
	 * @return  Success
	 */
	template<class T> bool
	Write (const Matrix<T>& M, const TiXmlElement* txe) {
		std::string uri (txe->Attribute("uri"));
		return this->Write (M, uri);
	}

private:

	NIFile ();
	NIFile (const NIFile&);


};

#define niwrite(X,Y) _niwrite (X,Y,#X)
	template<class T> inline static bool
	_niwrite (const Matrix<T>& M, const std::string& fname, const std::string& uri) {

		NIFile nif (fname, WRITE);
		nif.Write(M, uri);
		return true;

	}


	template<class T> inline static Matrix<T>
	niread (const std::string& fname, const std::string& uri = "") {

		NIFile nif (fname, READ);
		return nif.Read<T>(uri);

	}




#endif /* NIFILE_HPP_ */
