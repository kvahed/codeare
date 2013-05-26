/*
 * CODFile.hpp
 *
 *  Created on: May 26, 2013
 *      Author: kvahed
 */

#ifndef __CODFILE_HPP__
#define __CODFILE_HPP__

namespace codeare {
namespace matrix{
namespace io{

	class CODFile : public IOFile {



	public:


		CODFile (const std::string& fname, const IOMode mode,
				const Params& params = Params(), const bool verbose = false) :
					IOFile (fname, mode, params, verbose) {}


		~CODFile () {};


		template <class T> Matrix<T>
		Read (const std::string uri) {

		}

	private:


		CODFile () {};

		CODFile (const CODFile&) {};

	};

}
}
}



#endif /* __CODFILE_HPP__ */
