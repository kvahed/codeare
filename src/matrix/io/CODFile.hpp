/*
 * CODFile.hpp
 *
 *  Created on: May 26, 2013
 *      Author: kvahed
 */

#ifndef __CODFILE_HPP__
#define __CODFILE_HPP__


#include "Matrix.hpp"
#include "IOFile.hpp"

#include <fstream>
#include <iostream>


namespace codeare {
namespace matrix{
namespace io{

	template <class T>
	struct CODTraits;

	template<>
	struct CODTraits<float> {

	};

	static char* delimiter = "543f562189f1e82beb9c177f89f67822";

	class CODFile : public IOFile {

	public:


		CODFile (const std::string& fname, const IOMode mode = READ,
				const Params& params = Params(), const bool verbose = false) :
					IOFile (fname, mode, params, verbose) {



			const char* R = "rb";
			const char* W = "wb";

			bool  reading = (mode == READ);

			if (reading)
				assert (fexists(fname));

			if ((m_file = fopen(this->m_fname.c_str(), reading ? R : W))==NULL)
				printf("Cannot open %s file (%s).\n", this->m_fname.c_str(), reading ? R : W);

		}


		~CODFile () {
			if (m_file)
				fclose(m_file);
		};


		template <class T> Matrix<T>
		Read (const std::string& uri = "") const {

			int dt;
			size_t  ns, n;
			std::vector<size_t> dims (INVALID_DIM,1);
			float res[INVALID_DIM];
			char* name;
			Matrix<T> M;
			T t;

			// Read type
			if (!mread ( &dt, sizeof(   int),           1, m_file, "data type")) return M;

			// Matrix and data type must fit as of now.
			if ((typeid(T) == float_type    && dt == RLFL) ||
				(typeid(T) == double_type   && dt == RLDB) ||
				(typeid(T) == cxfl_type     && dt == CXFL) ||
				(typeid(T) == cxdb_type     && dt == CXDB) ||
				(typeid(T) == typeid(long)  && dt == LONG) ||
				(typeid(T) == typeid(short) && dt == SHRT)) {

				// Read dimensions and allocate matrix
				if (!mread (&dims[0], sizeof(size_t), INVALID_DIM, m_file, "dimensions")) return M;
				M = Matrix<T>(dims);
				n = numel(M);

				// Read resolutions and assign
				if (!mread ( res, sizeof( float), INVALID_DIM, m_file, "resolutions")) return M;
				for (size_t i = 0; i < INVALID_DIM; i++)
					M.Res(i) = res[i];

				// Name
				if (!mread ( &ns, sizeof(size_t),           1, m_file, "name length")) return M;
				name = new char [ns];
				if (!mread (name,   sizeof(char),          ns, m_file, "name")) return M;
				M.SetClassName(name);

				// Read data
				if (!mread (&M[0],     sizeof(T),           n, m_file, "data")) return M;

				//Close and clean up;
				delete name;

			}

			return M;

		}


		template <class T> bool
		Write (const Matrix<T>& M, const std::string& uri = "") {

			assert (m_file != NULL);

			size_t n = numel(M);
			int dt;
			size_t ns, l;

			if      (typeid(T) == typeid(float))  dt = RLFL;
			else if (typeid(T) == typeid(double)) dt = RLDB;
			else if (typeid(T) == typeid(cxfl))   dt = CXFL;
			else if (typeid(T) == typeid(cxdb))   dt = CXDB;
			else if (typeid(T) == typeid(long))   dt = LONG;
			else if (typeid(T) == typeid(short))  dt = SHRT;

			// Dump type
			if (!mwrite(&dt, sizeof(int), 1, m_file, "data type"))
				return false;

			// Dump dimensions
			if (!mwrite((const void*) M.Dim(), sizeof(size_t), INVALID_DIM, m_file, "dimensions"))
				return false;

			// Dump resolutions
			if (!mwrite((const void*) M.Res(), sizeof(float), INVALID_DIM, m_file, "resolutions"))
				return false;

			// Size of name and name
			ns = uri.size();
			if (!mwrite(&ns, sizeof(size_t), 1, m_file, "name length"))
				return false;

			// Dump name
			if (!mwrite((const void*) uri.c_str(), sizeof(char), ns, m_file, "name"))
				return false;

			// Dump data
			if (!mwrite((const void*) M.Memory(), sizeof(T), n, m_file, "data"))
				return false;

			return true;

		}

	private:


		CODFile () : m_file (0) {};
		CODFile (const CODFile&) : m_file(0) {};

		FILE* m_file;

	};


	template<class T>
	static bool codwrite (const std::string& fname, const Matrix<T>& M);

	template<class T>
	static Matrix<T> codread (const std::string& fname);

}
}
}



#endif /* __CODFILE_HPP__ */
