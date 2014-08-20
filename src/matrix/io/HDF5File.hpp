/*
 * HDF5File.hpp
 *
 *  Created on: Jan 17, 2013
 *      Author: kvahed
 */

#ifndef __HDF5FILE_HPP__
#define __HDF5FILE_HPP__

#include "IOFile.hpp"
#include "Matrix.hpp"
#include "Tokenizer.hpp"

#include <boost/tokenizer.hpp>

#include <H5Cpp.h>
using namespace H5;

namespace codeare {
namespace matrix {
namespace io {

	template<class T> struct HDF5Traits;

	template<> struct HDF5Traits<float> {
		static PredType* PType () {
			return (PredType*) new FloatType (PredType::NATIVE_FLOAT);
		}
	};
	template<> struct HDF5Traits<double> {
		static PredType* PType () {
			return (PredType*) new FloatType (PredType::NATIVE_DOUBLE);
		}
	};
	template<> struct HDF5Traits<cxfl> {
		static PredType* PType () {
			return (PredType*) new FloatType (PredType::NATIVE_FLOAT);
		}
	};
	template<> struct HDF5Traits<cxdb> {
		static PredType* PType () {
			return (PredType*) new FloatType (PredType::NATIVE_DOUBLE);
		}
	};
	template<> struct HDF5Traits<short> {
		static PredType* PType () {
			return (PredType*) new FloatType (PredType::NATIVE_SHORT);
		}
	};

	class HDF5File : public IOFile {

	public:

		/**
		 * @brief   Open HDF5 file
		 *
		 * @param  fname   File name
		 * @param  mode    IO mode (R/RW)
		 * @param  params  Optional params
		 * @param  verbose Verbosity
		 */
		HDF5File  (const std::string& fname, const IOMode mode = READ,
				Params params = Params(), const bool verbose = false) :
					IOFile(fname, mode, params, verbose) {

			Exception::dontPrint();

			try {
				m_file = H5File (fname, (mode == READ) ? H5F_ACC_RDONLY :H5F_ACC_TRUNC);
				if (this->m_verb)
					printf ("File %s opened %s\n", fname.c_str(), (mode == READ) ? "R" : "RW");
			} catch (const FileIException& e) {
				printf ("Opening %s failed\n", fname.c_str());
				e.printError();
			}

			this->m_status = OK;

		}



		/**
		 * @brief  Default destructor
		 */
		virtual ~HDF5File () {
			Close ();
		}



		/**
		 * @brief   Clean up and close file
		 */
		virtual void
		Close () {
			try {
				m_file.flush(H5F_SCOPE_LOCAL);
			} catch (const Exception& e) {
				this->m_status = HDF5_ERROR_FFLUSH;
				printf ("Couldn't flush HDF5 file %s!\n%s\n", this->FileName().c_str(), e.getDetailMsg().c_str());
			}
			try {
				m_file.close();
			} catch (const Exception& e) {
				this->m_status = HDF5_ERROR_FCLOSE;
				printf ("Couldn't close HDF5 file %s!\n%s\n", this->FileName().c_str(), e.getDetailMsg().c_str());
			}
		}



		template<class T> Matrix<T>
		Read (const std::string& uri) const throw () {

			T         t       = (T) 0;
			DataSet   dataset = m_file.openDataSet(uri);
			DataSpace space   = dataset.getSpace();
			Vector<hsize_t> dims (space.getSimpleExtentNdims());
			size_t    ndim    = space.getSimpleExtentDims(&dims[0], NULL);

			if (this->m_verb) {
				printf ("Reading dataset %s ... ", uri.c_str());
				fflush(stdout);
			}

			if (is_complex(t)) {
				dims.PopBack();
				--ndim;
			}

			Vector<size_t> mdims (ndim,1);

			for (size_t i = 0; i < ndim; ++i)
				mdims[i] = dims[ndim-i-1];

			PredType* type = HDF5Traits<T>::PType();

			Matrix<T> M (mdims);
			dataset.read (&M[0], *type);

			if (this->m_verb)
				printf ("O(%s) done\n", DimsToCString(M));

			space.close();
			dataset.close();

			return M;

		}



		template<class T> bool
		Write (const Matrix<T>& M, const std::string& uri) throw () {

			T t = (T)0;
			Group group, *tmp;
			std::string path;

            boost::tokenizer<> tok(uri);
			std::vector<std::string> sv (Split (uri, "/"));
			std::string name = sv[sv.size() - 1];
			sv.pop_back(); // data name not part of path

			if (sv.size() == 0)
				path = "/";
			else
				for (size_t i = 0; i < sv.size(); i++) {
					if (sv[i].compare(""))
						path += "/";
						path += sv[i];
				}

			if (this->m_verb)
				printf ("Creating dataset %s at path (%s)\n", name.c_str(), path.c_str());

			try {

				group = m_file.openGroup(path);
				if (this->m_verb)
					printf ("Group %s opened for writing\n", path.c_str()) ;

			} catch (const Exception&) {

				for (size_t i = 0, depth = 0; i < sv.size(); i++) {

					if (sv[i].compare("")) {

						try {
							group = (depth) ? (*tmp).openGroup(sv[i])   : m_file.openGroup(sv[i]);
						} catch (const Exception&) {
							group = (depth) ? (*tmp).createGroup(sv[i]) : m_file.createGroup(sv[i]);
						}

						tmp = &group;
						depth++;

					}

				}

			}

			// One more field for complex numbers
			size_t tmpdim = ndims(M);

			Vector<hsize_t> dims (tmpdim);

			for (size_t i = 0; i < tmpdim; i++)
				dims[i] = M.Dim(tmpdim-1-i);

			if (is_complex(t)) {
				dims.push_back(2);
				tmpdim++;
			}

			DataSpace space (tmpdim, &dims[0]);
			PredType*  type = HDF5Traits<T>::PType();

			DataSet set = group.createDataSet(name, (*type), space);

			set.write   (M.Ptr(), (*type));
			set.close   ();
			space.close ();

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

		HDF5File (const HDF5File&) {}
		HDF5File  () {}

		H5File m_file; /// @brief My file


	};

}// namespace io
}// namespace matrix
}// namespace codeare


	template<class T> inline static bool
	_h5write (const Matrix<T>& M, const std::string& fname, const std::string& uri) {
		using namespace codeare::matrix::io;
		HDF5File h5f (fname, WRITE);
		h5f.Write(M, uri);
		return true;
	}
#define h5write(X,Y) _h5write (X,Y,#X)


	template<class T> inline static Matrix<T>
	h5read (const std::string& fname, const std::string& uri) {
		using namespace codeare::matrix::io;
		HDF5File h5f (fname, READ);
		return h5f.Read<T>(uri);
	}





#endif /* __HDF5FILE_HPP__ */
