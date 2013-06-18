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

	template<> struct CODTraits<float> {
		static const dtype dt = RLFL;
	};
	template<> struct CODTraits<double> {
		static const dtype dt = RLDB;
	};
	template<> struct CODTraits<cxfl> {
		static const dtype dt = CXFL;
	};
	template<> struct CODTraits<cxdb> {
		static const dtype dt = CXDB;
	};
	template<> struct CODTraits<long> {
		static const dtype dt = LONG;
	};
	template<> struct CODTraits<short> {
		static const dtype dt = SHRT;
	};

	static const std::string delim = "543f562189f1e82beb9c177f89f67822";

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
			std::vector<size_t> dim;
			std::vector<float>  res;
			char* name;
			Matrix<T> M;
			T t;

			// Read type
			if (!mread ( &dt,           1, m_file, "data type")) return M;

			// Matrix and data type must fit as of now.
			if (CODTraits<T>::dt == dt) {

				if (!mread (&n, 1, m_file, "dimensions")) return M;
				dim = std::vector<size_t>(n,1);
				res = std::vector<float>(n,1.0);

				// Read dimensions and allocate matrix
				if (!mread (dim, m_file, "dimensions")) return M;

				//Read resolutions and assign
				if (!mread (res, m_file, "resolutions")) return M;
				M = Matrix<T>(dim,res);
				n = numel(M);

				// Name
				if (!mread (&n,  1, m_file, "name length")) return M;

				name = new char [n+1];

				if (!mread (name, n, m_file, "name")) return M;
				name[n] = '\0';
				M.SetClassName(name);

				// Read data
				if (!mread (M.Container(), m_file, "data")) return M;

				//Close and clean up;
				delete name;

			}

			return M;

		}


		template <class T> bool
		Write (const Matrix<T>& M, const std::string& uri = "") {

			assert (m_file != NULL);

			dtype dt = CODTraits<T>::dt;
			size_t n, l;

			// Dump type
			if (!mwrite(&dt,         1, m_file, "data type"))
				return false;

			n = M.NDim();
			if (!mwrite(&n,          1, m_file, "data dimensions"))
				return false;

			// Dump dimensions
			if (!mwrite(M.Dim(),        m_file, "dimensions"))
				return false;

			// Dump resolutions
			if (!mwrite(M.Res(),        m_file, "resolutions"))
				return false;

			// Size of name and name
			n = uri.size();
			if (!mwrite(&n,          1, m_file, "name length"))
				return false;

			// Dump name
			if (!mwrite(uri.c_str(), n, m_file, "name"))
				return false;

			// Dump data
			if (!mwrite(M.Container(),  m_file, "data"))
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
