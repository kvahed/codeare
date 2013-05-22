/*
 * MATFile.hpp
 *
 *  Created on: May 20, 2013
 *      Author: kvahed
 */

#ifndef MATFILE_HPP_
#define MATFILE_HPP_

#include "IOFile.hpp"
#include "mat.h"

#ifdef HAVE_CXXABI_H
    #include <cxxabi.h>
#endif


static inline std::string
deman (const char* symbol) {

#ifdef HAVE_CXXABI_H

	size_t size;
	int    status;
	char   temp[128];
	char*  demangled;

	//first, try to demangle a c++ name
	if (1 == sscanf(symbol, "%*[^(]%*[^_]%127[^)+]", temp)) {

		if (NULL != (demangled = abi::__cxa_demangle(temp, NULL, &size, &status))) {
			std::string result(demangled);
			free(demangled);
			return result;
		}
	}

	//if that didn't work, try to get a regular c symbol
	if (1 == sscanf(symbol, "%127s", temp)) {
		return temp;
	}

#endif

	//if all else fails, just return the symbol
	return symbol;

}

template <class T> inline static bool
MXValidate  (const Matrix<T>& M, const mxArray* mxa) {

	T t;
	mxClassID     mcid = mxGetClassID(mxa);
	std::string cplx = (mxIsComplex(mxa)) ? "complex" : "real";

	const char* vname = deman(typeid(T).name()).c_str();

	if (is_singlep(t) && mcid == 7)
		return true;
	if (is_doublep(t) && mcid == 6)
		return true;
	else {
		printf ("Matrix is %s, yet Matlab variable is %s %s\n", vname, mxGetClassName(mxa), cplx.c_str());
		return false;
	}

	return true;

}

namespace codeare {
namespace matrix {
namespace io {

	class MLFile : public IOFile {

	public:

		MLFile (const std::string& fname, const IOMode mode, const Params& params = Params(), const bool verbose = false) :
			IOFile (fname, mode, params, verbose) {

			if (verbose)
				printf ("Opening %s\n", fname.c_str());
			m_file = matOpen (fname.c_str(), (mode == READ) ? "r" : "w" );

			if (m_file == NULL)
				printf ("Error opening MATLAB file %s\n", fname.c_str());

		}

		virtual ~MLFile () {

			if (m_file)
				if (matClose(m_file) != 0)
				printf ("Error closing file %s\n",this->m_fname.c_str());

		}

		template<class T> Matrix<T>
		Read (const std::string& uri) const {

			mxArray*      mxa = matGetVariable (m_file, uri.c_str());
			Matrix<T> M;

			if (!mxa) {
				printf ("Error opening variable %s\n", uri.c_str());
				return M;
			}

			mxClassID     mcid = mxGetClassID(mxa);
			mwSize        ndim = mxGetNumberOfDimensions(mxa);
			const mwSize*  dim = mxGetDimensions(mxa);
			size_t i = 0;


			if (!MXValidate (M, mxa))
				return M;

			std::vector<size_t> mdims(ndim,1);

			for (; i < ndim; i++)
				mdims[i] = (size_t)dim[i];

			M = Matrix<T> (mdims);
			// -------------------------------------------

			// Copy from memory block ----------------------

			printf ("  Reading (type: %i) %s: (%s) ... ", (int)mcid, uri.c_str(), DimsToCString(M)); fflush(stdout);

			if        (typeid(T) == typeid(cxfl)) {
				if (mxIsComplex(mxa))
					for (i = 0; i < numel(M); i++) {
						float f[2] = {((float*)mxGetPr(mxa))[i], ((float*)mxGetPi(mxa))[i]}; // Template compilation. Can't create T(real,imag)
						memcpy(&M[i], f, 2 * sizeof(float));
					}
				else
					for (i = 0; i <numel(M); i++) {
						float f[2] = {((float*)mxGetPr(mxa))[i], 0.0};
						memcpy(&M[i], f, 2 * sizeof(float));
					}
			} else if (typeid(T) == typeid(cxdb)) {
				if (mxIsComplex(mxa))
					for (i = 0; i < numel(M); i++) {
						double* tmp = (double*) &M[0];
						tmp [i*2]   = mxGetPr(mxa)[i];
						tmp [i*2+1] = mxGetPi(mxa)[i];
					}
				else
					for (i = 0; i <numel(M); i++) {
						double* tmp = (double*) &M[0];
						tmp [i*2]   = mxGetPr(mxa)[i];
					}
			} else
				memcpy (&M[0], mxGetPr(mxa), numel(M) * sizeof (T));


			printf ("done\n");
			// -------------------------------------------

			// Clean up and close file -------------------
			if (mxa != NULL)
				mxDestroyArray(mxa);

			// -------------------------------------------

			return M;


		}

		template <class T> bool
		Write (const Matrix<T>& M, const std::string& uri) {


			// Declare dimensions and allocate array -----

			mwSize   dim[INVALID_DIM];

			for (size_t i = 0; i < INVALID_DIM; i++)
				dim[i] = (mwSize)M.Dim(i);

			mxArray*  mxa = 0;

			if      (typeid(T) == float_type)
				mxa = mxCreateNumericArray (INVALID_DIM, dim, mxSINGLE_CLASS,    mxREAL);
			else if (typeid(T) == double_type)
				mxa = mxCreateNumericArray (INVALID_DIM, dim, mxDOUBLE_CLASS,    mxREAL);
			else if (typeid(T) == cxfl_type)
				mxa = mxCreateNumericArray (INVALID_DIM, dim, mxSINGLE_CLASS, mxCOMPLEX);
			else if (typeid(T) == cxdb_type)
				mxa = mxCreateNumericArray (INVALID_DIM, dim, mxDOUBLE_CLASS, mxCOMPLEX);
			else if (typeid(T) == typeid(short))
				mxa = mxCreateNumericArray (INVALID_DIM, dim,  mxINT16_CLASS,    mxREAL);
			// -------------------------------------------


			// Copy to memory block ----------------------

			if (typeid(T) == cxfl_type) {
			    float* re = (float*)mxGetData(mxa);
				float* im = (float*)mxGetImagData(mxa);
				for (size_t i = 0; i < numel(M); i++) {
					re[i] = creal(M[i]);
					im[i] = cimag(M[i]);
				}
			} else if (typeid(T) == typeid(cxdb)) {
				double* re = mxGetPr(mxa);
				double* im = mxGetPi(mxa);
				for (size_t i = 0; i < numel(M); i++) {
					re[i] = creal(M[i]);
					im[i] = cimag(M[i]);
				}
			} else
				memcpy(mxGetData(mxa), M.Memory(), numel(M) * sizeof(T));

			// -------------------------------------------

			// Write data --------------------------------
			int status = matPutVariable(m_file, uri.c_str(), mxa);

			if (status != 0) {
				printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
				return false;
			}
			// -------------------------------------------


			// Clean up RAM ------------------------------

			if (mxa != NULL)
				mxDestroyArray(mxa);
			// -------------------------------------------

			return true;


		}


		/**
		 * @brief Read a particular data set from file
		 *
		 * @return  Success
		 */
		template<class T> Matrix<T>
		Read (const TiXmlElement* txe) const {
			std::string uri (txe->Attribute("uri"));
			return this->Read<T>(uri);
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

		MATFile*  m_file;

	};

}}}


#endif /* MATFILE_HPP_ */
