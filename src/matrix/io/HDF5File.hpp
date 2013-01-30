/*
 * HDF5File.hpp
 *
 *  Created on: Jan 17, 2013
 *      Author: kvahed
 */

#ifndef __HDF5FILE_HPP__
#define __HDF5FILE_HPP__

#include "IOFile.hpp"

namespace codeare {
namespace matrix {
namespace io {

	class HDF5File : public IOFile {

		HDF5File  (const std::string& fname, const IOMode mode, const bool verbosity) {}

		virtual ~HDF5File () {}

	private:

		HDF5File (const HDF5File& h5f) {}

		HDF5File  () {}

	};


}// namespace io
}// namespace matrix
}// namespace codeare


#endif /* __HDF5FILE_HPP__ */
