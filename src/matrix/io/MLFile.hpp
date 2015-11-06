/*
 * MATFile.hpp
 *
 *  Created on: May 20, 2013
 *      Author: kvahed
 */

#ifndef MATFILE_HPP_
#define MATFILE_HPP_

#include "Workspace.hpp"
#include "IOFile.hpp"
#include "mat.h"

#include "Demangle.hpp"

template <class T> inline static bool MXValidate  (const Matrix<T>& M, const mxArray* mxa) {

	T t = (T)0;
	mxClassID     mcid = mxGetClassID(mxa);
	std::string cplx = (mxIsComplex(mxa)) ? "complex" : "real";

	const char* vname = demangle(typeid(T).name()).c_str();

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

template <class T>
struct MXTraits;

template <class T> static void write_real (mxArray* mxa, const Matrix<T>& M) {
	memcpy(mxGetData(mxa), M.Ptr(), numel(M) * sizeof(T));
}

template <class T> static void read_real (Matrix<T>& M, const mxArray* mxa) {
	memcpy (&M[0], mxGetPr(mxa), numel(M) * sizeof (T));
}

template <class T> static void write_complex (mxArray* mxa, const Matrix<std::complex<T> >& M) {
	T* re = (T*)mxGetPr(mxa);
	T* im = (T*)mxGetPi(mxa);
	for (size_t i = 0; i < numel(M); ++i) {
		re[i] = real(M[i]);
		im[i] = imag(M[i]);
	}
}

template <class T> static void read_complex (Matrix<std::complex<T> >& M, const mxArray* mxa) {
	T* re = (T*)mxGetPr(mxa);
	T* im = (T*)mxGetPi(mxa);
	bool cplx = (im != NULL);
	for (size_t i = 0; i < numel(M); ++i)
		M[i] = std::complex<T>(re[i], (cplx) ? im[i] : 0.0);
}

template <> struct MXTraits<float> {
	static const mxClassID prec = mxSINGLE_CLASS;
	static const mxComplexity cplx = mxREAL;
	typedef float T;
	static void Write (mxArray* mxa, const Matrix<float>& M) {write_real(mxa,M);}
	static void Read (Matrix<T>& M, const mxArray* mxa) {read_real(M,mxa);}
};
template <> struct MXTraits<double> {
	static const mxClassID prec = mxDOUBLE_CLASS;
	static const mxComplexity cplx = mxREAL;
	typedef double T;
	static void Write (mxArray* mxa, const Matrix<T>& M) {write_real(mxa,M);}
	static void Read (Matrix<T>& M, const mxArray* mxa) {read_real(M,mxa);}
};
template <> struct MXTraits<cxfl> {
	static const mxClassID prec = mxSINGLE_CLASS;
	static const mxComplexity cplx = mxCOMPLEX;
	typedef cxfl T;
	typedef float T2;
	static void Write (mxArray* mxa, const Matrix<T>& M) {write_complex(mxa,M);}
	static void Read (Matrix<T>& M, const mxArray* mxa) {read_complex(M,mxa);}
};
template <> struct MXTraits<cxdb> {
	static const mxClassID prec = mxDOUBLE_CLASS;
	static const mxComplexity cplx = mxCOMPLEX;
	typedef cxdb T;
	typedef double T2;
	static void Write (mxArray* mxa, const Matrix<T>& M) {write_complex(mxa,M);}
	static void Read (Matrix<T>& M, const mxArray* mxa) {read_complex(M,mxa);}
};
template <> struct MXTraits<short> {
	static const mxClassID prec = mxINT16_CLASS;
	static const mxComplexity cplx = mxREAL;
	typedef cxdb T;
	static void Write (mxArray* mxa, const Matrix<T>& M) {write_complex(mxa,M);}
	static void Read (Matrix<T>& M, const mxArray* mxa) {read_real(M,mxa);}
};

	class MLFile : public IOFile {

	public:

		MLFile (const std::string& fname, const IOMode mode, const Params& params = Params(),
				const bool verbose = false) : IOFile (fname, mode, params, verbose) {
			if (verbose)
				printf ("Opening %s (%s).\n", fname.c_str(), (mode == READ) ? "r" : "w");
			m_file = matOpen (fname.c_str(), (mode == READ) ? "r" : "w" );
			if (m_file == NULL) {
				printf ("Error opening MATLAB file %s\n", fname.c_str());
				assert (false);
			}
		}

		virtual ~MLFile () {
			if (m_file)
				if (matClose(m_file) != 0)
				printf ("Error closing file %s\n",this->m_fname.c_str());
		}

		template<class T> Matrix<T>	Read (const std::string& uri) const {

			mxArray*      mxa = matGetVariable (m_file, uri.c_str());
			Matrix<T> M;

			if (!mxa) {
				printf ("**ERROR**: Failed to retrieve variable %s\n", uri.c_str());
				assert (false);
			}

			mxClassID     mcid = mxGetClassID(mxa);
			mwSize        ndim = mxGetNumberOfDimensions(mxa);
			const mwSize*  dim = mxGetDimensions(mxa);
			size_t i = 0;

			assert (MXValidate (M, mxa));

			Vector<size_t> mdims(ndim,1);

			for (; i < ndim; ++i)
				mdims[i] = (size_t)dim[i];

			M = Matrix<T> (mdims);

			// Copy from memory block ----------------------
			MXTraits<T>::Read(M, mxa);

			// Clean up and close file -------------------
			if (mxa != NULL)
				mxDestroyArray(mxa);

			return M;


		}

		template <class T> bool	Write (const Matrix<T>& M, const std::string& uri) {

			// Declare dimensions and allocate array
			size_t nd = M.NDim();
			Vector<mwSize> dim(nd);
			for (size_t i = 0; i < nd; ++i)
				dim[i] = (mwSize)M.Dim(i);
			mxArray*  mxa = mxCreateNumericArray (nd, &dim[0], MXTraits<T>::prec, MXTraits<T>::cplx);

			// Assign data
			MXTraits<T>::Write (mxa, M);

			// Write data
			int status = matPutVariable(m_file, uri.c_str(), mxa);
			if (status != 0) {
				printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
				return false;
			}

			// Clean up RAM
			if (mxa != NULL)
				mxDestroyArray(mxa);

			return true;


		}


		/**
		 * @brief Read a particular data set from file
		 *
		 * @return  Success
		 */
		template<class T> Matrix<T>	Read (const TiXmlElement* txe) const {
			std::string uri (txe->Attribute("uri"));
			return this->Read<T>(uri);
		}


		/**
		 * @brief  Write data to file
		 *
		 * @return  Success
		 */
		template<class T> bool Write (const Matrix<T>& M, const TiXmlElement* txe) {
			std::string uri (txe->Attribute("uri"));
			return this->Write (M, uri);
		}


        inline void Read () const {
            std::cout  << "    File name: " << m_fname << std::endl;
        	int vars = 0;
        	char** dir = matGetDir(m_file, &vars);
        	for (auto i = 0; i < vars; ++i)
        		DoDataset (dir[i]);
        }

        inline void DoDataset (char* name) const {
        	std::cout << "      Dataset: " << name << ": ("; fflush(stdout);
    		mxArray* mxa = matGetVariable (m_file, name);
			mwSize        ndim = mxGetNumberOfDimensions(mxa);
			const mwSize*  dim = mxGetDimensions(mxa);
			for (mwSize i=0; i < ndim-1; ++i)
				std::cout << dim[i] << " ";
			std::cout << dim[ndim-1] << ") [";
			switch (mxGetClassID(mxa))
			{
				case mxUNKNOWN_CLASS:
					std::cout << "unknown";
					break;
				case mxSINGLE_CLASS:
					if (mxIsComplex(mxa)) {
						std::cout << "complex ";
						Matrix<cxfl> M = Read<cxfl>(name);
						wspace.Add(name, M);
					} else {
						Matrix<float> M = Read<float>(name);
						wspace.Add(name, M);
					}
					std::cout << "float";
					break;
				case mxDOUBLE_CLASS:
					if (mxIsComplex(mxa)) {
						std::cout << "complex ";
						Matrix<cxdb> M = Read<cxdb>(name);
						wspace.Add(name, M);
					} else {
						Matrix<double> M = Read<double>(name);
						wspace.Add(name, M);
					}
					std::cout << "double";
					break;
		        case mxCELL_CLASS:
					std::cout << "cell";
					break;
				case mxSTRUCT_CLASS:
					std::cout << "struct";
					break;
				case mxLOGICAL_CLASS:
					std::cout << "logical";
					break;
				case mxCHAR_CLASS:
					std::cout << "char";
					break;
				case mxVOID_CLASS:
					std::cout << "void";
					break;
				case mxINT8_CLASS:
					std::cout << "int8_t";
					break;
				case mxUINT8_CLASS:
					std::cout << "uint8_t";
					break;
				case mxINT16_CLASS:
					std::cout << "int16_t";
					break;
				case mxUINT16_CLASS:
					std::cout << "uint16_t";
					break;
				case mxINT32_CLASS:
					std::cout << "int32_t";
					break;
				case mxUINT32_CLASS:
					std::cout << "uint32_t";
					break;
				case mxINT64_CLASS:
					std::cout << "int64_t";
					break;
				case mxUINT64_CLASS:
					std::cout << "uint16_t";
					break;
				case mxFUNCTION_CLASS:
					std::cout << "function";
					break;
				default:
					std::cout << "unknown";
					break;
			}

			std::cout << "]" << std::endl;

        }

	private:

		MLFile (const MLFile& mf) : m_file(0) {}
		MLFile ()                 : m_file(0) {}
		MATFile*                    m_file;

	};


#define mxwrite(X,Y) _mxwrite (X,Y,#X)
	template<class T> inline static bool
	_mxwrite (const Matrix<T>& M, const std::string& fname, const std::string& uri) {

		MLFile mf (fname, WRITE);
		mf.Write(M, uri);
		return true;

	}


	template<class T> inline static Matrix<T>
	mxread (const std::string& fname, const std::string& uri) {

		MLFile mf (fname, READ);
		return mf.Read<T>(uri);

	}

}}}


#endif /* MATFILE_HPP_ */
