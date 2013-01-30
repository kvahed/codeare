/*
 * IOContext.hpp
 *
 *  Created on: Jan 16, 2013
 *      Author: kvahed
 */

#ifndef __IOCONTEXT_HPP__
#define __IOCONTEXT_HPP__

#include "ISMRMRD.hpp"
#include "HDF5File.hpp"
#include "SyngoFile.hpp"
/*
#include "NetCDF.hpp"
#include "Nifti.hpp"
#include "CodRaw.hpp"
#include "GE.hpp"
#include "PHILIPS.hpp"
*/
namespace codeare {
namespace matrix{
namespace io{

	/**
	 * @brief Supported data formats
	 */
	enum IOStrategy {CODRAW, HDF5, MATLAB, ISMRM, SYNGOMR, GE, PHILIPS};

	/**
	 * @brief       Interface to concrete IO implementation
	 */
	class IOContext {

	public:

		/**
		 * @brief   Default constructor
		 */
		IOContext () :
			m_iof(0) {}

		/**
		 * @brief   Default copy constructor
		 */
		IOContext (const IOContext& ioc) :
			m_iof (0) {}

		/**
		 *
		 */
		IOContext (const std::string& fname, const IOStrategy& iocs = HDF5,
				const Params& params = Params(), const IOMode& mode = READ,
				const bool& verbosity = true) :
			m_iof(0) {

			switch (iocs) {
				case CODRAW:  break;
				case HDF5:    break;
				case MATLAB:  break;
				case ISMRM: m_iof = (IOFile*) new IRDFile (fname, mode, params, verbosity);break;
				case SYNGOMR: m_iof = (IOFile*) new SyngoFile (fname, mode, params, verbosity); break;
				case GE:      break;
				case PHILIPS: break;
				default:      break;
			}

		}

		/**
		 * @brief  Destroy file handle
		 */
		~IOContext () {
			if (m_iof)
				delete m_iof;
		}

		IOFile* Context () {
			return m_iof;
		}

		/**
		 * @brief   Return concrete handle's status
		 */
		RRSModule::error_code Status () const {

			if (m_iof)
				return m_iof->Status();

			return RRSModule::NULL_FILE_HANDLE;

		}

		/**
		 * @brief   Return concrete handle's status
		 */
		RRSModule::error_code Read (const std::string& dname) const {

			if (m_iof)
				return m_iof->Read(dname);

			return RRSModule::NULL_FILE_HANDLE;

		}

		/**
		 * @brief   Return concrete handle's status
		 */
		RRSModule::error_code Write (const std::string& dname) const {

			if (m_iof)
				return m_iof->Write(dname);

			return RRSModule::NULL_FILE_HANDLE;

		}


		/**
		 * @brief   Return concrete handle's status
		 */
		RRSModule::error_code Read (const char* dname) const {

			if (m_iof)
				return m_iof->Read(std::string(dname));

			return RRSModule::NULL_FILE_HANDLE;

		}

		/**
		 * @brief   Return concrete handle's status
		 */
		RRSModule::error_code Write (const char* dname) const {

			if (m_iof)
				return m_iof->Write(std::string(dname));

			return RRSModule::NULL_FILE_HANDLE;

		}

		/**
		 * @brief   Return concrete handle's status
		 */
		boost::any& GetAny (const std::string& dname) const {

			if (!m_iof)
				return Toolbox::Instance()->void_any;

			return m_iof->GetAny (dname);

		}


		/**
		 * @brief   Return concrete handle's status
		 */
		template<class T>
		const Matrix<T>& Get (const std::string& dname) const {

			return m_iof->Get<T> (dname);
		}


		/**
		 * @brief   Return concrete handle's status
		 */
		RRSModule::error_code CleanUp () const {

			if (m_iof)
				return m_iof->CleanUp();

			return RRSModule::NULL_FILE_HANDLE;

		}


		IOFile* m_iof; /**< @brief  My actual context */

	};

}
}
}


#endif /* __IOCONTEXT_HPP__ */
