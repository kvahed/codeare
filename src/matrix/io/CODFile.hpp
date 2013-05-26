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

	class CODFile : public IOFile {



	public:


		CODFile (const std::string& fname, const IOMode mode = READ,
				const Params& params = Params(), const bool verbose = false) :
					IOFile (fname, mode, params, verbose) {

			const char* R = "rb";
			const char* W = "wb";

			assert (fexists(fname));
			if ((m_file = fopen(this->m_fname.c_str(), (mode == READ) ? R : W))==NULL)
				printf("Cannot open %s file (%s).\n", this->m_fname.c_str(), (mode == READ) ? R : W);

		}


		~CODFile () {
			if (m_file)
				fclose(m_file);
		};


		template <class T> Matrix<T>
		Read (const std::string& uri = "") const {

			Matrix<T> M;
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
	static bool codwrite (const std::string& fname, const Matrix<T>& M)

}
}
}



#endif /* __CODFILE_HPP__ */
