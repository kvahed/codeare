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

#include <H5Cpp.h>
using namespace H5;

namespace codeare {
namespace matrix {
namespace io {


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
				printf ("Couldn't flush HDF5 file %s!\n", this->FileName().c_str());
			}
			try {
				m_file.close();
			} catch (const Exception& e) {
				this->m_status = HDF5_ERROR_FCLOSE;
				printf ("Couldn't close HDF5 file %s!\n", this->FileName().c_str());
			}
		}



		template<class T> Matrix<T>
		Read (const std::string& uri) const {

			T         t;
			DataSet   dataset = m_file.openDataSet(uri);
			DataSpace space   = dataset.getSpace();
			hsize_t*  dims    = (hsize_t*) malloc (space.getSimpleExtentNdims() * sizeof (hsize_t));
			size_t    ndim    = MIN(space.getSimpleExtentDims(dims, NULL), (is_complex(t)) ? INVALID_DIM + 1 : INVALID_DIM);

			if (this->m_verb)
				printf ("Reading dataset %s ... ", uri.c_str());
			fflush(stdout);

			if (is_complex(t))
				--ndim;

			std::vector<size_t> mdims (ndim,1);

			for (size_t i = 0; i < ndim; ++i)
				mdims[i] = dims[ndim-i-1];

			PredType* type;
			if      (is_singlep(t))
				type = (PredType*) new FloatType (PredType::NATIVE_FLOAT);
			else if (is_doublep(t))
				type = (PredType*) new FloatType (PredType::NATIVE_DOUBLE);
			else
				type = (PredType*) new IntType   (PredType::NATIVE_SHORT);

			Matrix<T> M (mdims);
			dataset.read (&M[0], *type);


			if (this->m_verb)
				printf ("O(%s) done\n", DimsToCString(M));

			space.close();
			dataset.close();

			return M;

		}



		template<class T> bool
		Write (const Matrix<T>& M, const std::string& uri) {

			T t;
			Group group, *tmp;
			std::string path;

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

			} catch (const Exception& e) {

				for (size_t i = 0, depth = 0; i < sv.size(); i++) {

					if (sv[i].compare("")) {

						try {
							group = (depth) ? (*tmp).openGroup(sv[i])   : m_file.openGroup(sv[i]);
						} catch (const Exception& e) {
							group = (depth) ? (*tmp).createGroup(sv[i]) : m_file.createGroup(sv[i]);
						}

						tmp = &group;
						depth++;

					}

				}

			}

			// One more field for complex numbers
			size_t tmpdim = (is_complex(t)) ? INVALID_DIM+1 : INVALID_DIM;
			hsize_t* dims = new hsize_t[tmpdim];

			for (size_t i = 0; i < INVALID_DIM; i++)
				dims[i] = M.Dim(INVALID_DIM-1-i);

			if (typeid(T) == cxfl_type || typeid(T) == cxdb_type)
				dims[INVALID_DIM] = 2;

			DataSpace space (tmpdim, dims);
			PredType*  type;

			delete [] dims;

			if      (is_singlep(t))
				type = (PredType*) new FloatType (PredType::NATIVE_FLOAT);
			else if (is_doublep(t))
				type = (PredType*) new FloatType (PredType::NATIVE_DOUBLE);
			else if (typeid(T) == typeid(short))
				type = (PredType*) new IntType   (PredType::NATIVE_SHORT);

			DataSet set = group.createDataSet(name, (*type), space);

			set.write   (M.Memory(), (*type));
			set.close   ();
			space.close ();

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

		HDF5File (const HDF5File& h5f) {}

		HDF5File  () {}

		H5File m_file; /// @brief My file


	};


	template<class T> inline static
	bool h5write (const Matrix<T>& M, const std::string& fname, const std::string& uri) {

		HDF5File file (fname, WRITE);
		file.Write(M, uri);

		return true;

	}

	template<class T> inline static
	bool h5read (const Matrix<T>& M, const std::string& fname, const std::string& uri) {

		HDF5File h5f (fname, READ);
		M = h5f.Read<T>(uri);

		return true;

	}



}// namespace io
}// namespace matrix
}// namespace codeare


#endif /* __HDF5FILE_HPP__ */
