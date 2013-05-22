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
#include "MLFile.hpp"

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
	enum IOStrategy {CODRAW = 0, HDF5, MATLAB, ISMRM, NIFTI, SYNGOMR, GE, PHILIPS, NO_STRATEGY};
	static const std::string IOSName[(int)NO_STRATEGY] = {"CODRAW", "HDF5", "MATLAB", "NIFTI", "SYNGOMR", "GE", "PHILIPS"};

	/**
	 * @brief       Interface to concrete IO implementation
	 */
	class IOContext {

	public:

		/**
		 * @brief   Default constructor
		 */
		IOContext () :
			m_iof(0), m_ios (HDF5) {}

		/**
		 * @brief   Default copy constructor
		 */
		IOContext (const IOContext& ioc) :
			m_iof (0), m_ios (HDF5) {}


		/**
		 *
		 */
		IOContext (const std::string& fname, const std::string& fmode) :
			m_iof(0), m_ios(HDF5) {

			assert (fmode.compare("r") == 0 || fmode.compare("rb") == 0 || fmode.compare("w") == 0 || fmode.compare("wb") == 0);

			if (HasSuffix (fname, ".h5"))
				m_ios = HDF5;
			else if (HasSuffix (fname, ".mat"))
				m_ios = MATLAB;
			else if (HasSuffix (fname, ".cod"))
				m_ios = CODRAW;
			else if (HasSuffix (fname, ".nii"))
				m_ios = NIFTI;
			else if (HasSuffix (fname, ".dat"))
				m_ios = SYNGOMR;

			IOMode mode;
			if (fmode.compare("r") == 0 || fmode.compare("rb") == 0)
				mode = READ;
			else if (fmode.compare("w") == 0 || fmode.compare("wb") == 0)
				mode = WRITE;

			this->Concretize(fname, mode, Params(), false);

		}



		/**
		 *
		 */
		IOContext (const std::string& fname, const IOStrategy& ios = HDF5,
				const IOMode& mode = READ, const Params& params = Params(),
				const bool verbosity = true) :
			m_iof(0) , m_ios(ios) {

			this->Concretize(fname, mode, params, verbosity);

		}


		/**
		 *
		 */
		IOContext (const TiXmlElement* txe, const std::string& base = ".",
				const IOMode mode = READ, const bool verbosity = true) :
			m_iof(0), m_ios(HDF5) {

			if (!txe)
				printf ("Ouch, XML element for construction is 0!\n");

			std::string fname = base + "/" + std::string(txe->Attribute("fname"));
			std::string ftype (txe->Attribute("ftype"));

			if (ftype.compare("CODRAW") == 0)
				m_ios = CODRAW;
			else if (ftype.compare("MATLAB") == 0)
				m_ios = MATLAB;
			else if (ftype.compare("HDF5") == 0)
				m_ios = HDF5;
			else if (ftype.compare("SYNGOMR") == 0)
				m_ios = SYNGOMR;

			this->Concretize(fname, mode, Params(), verbosity);

		} // TODO: UPPERCASE


		/**
		 * @brief  Destroy file handle
		 */
		~IOContext () {
			if (m_iof)
				delete m_iof;
		}

		IOFile& Context () {
			return *m_iof;
		}

		/**
		 * @brief   Return concrete handle's status
		 */
		error_code Status () const {

			if (m_iof)
				return m_iof->Status();

			return GENERAL_IO_ERROR;

		}

		/**
		 * @brief   Return concrete handle's status
		 */
		template<class T> Matrix<T>
		Read (const std::string& uri) const {

			if (m_iof) {
				switch (m_ios) 
					{
					case CODRAW:  break;
					case HDF5:    return ((HDF5File*)m_iof)->Read<T>(uri);
					case MATLAB:  return ((MLFile*)m_iof)->Read<T>(uri);
					case ISMRM:   return ((IRDFile*)m_iof)->Read<T>(uri);
					case SYNGOMR: return ((SyngoFile*)m_iof)->Read<T>(uri);
					case GE:      break;
					case PHILIPS: break;
					default:      break;
				}
			}

		}

		/**
		 * @brief   Return concrete handle's status
		 */
		template <class T>
		bool Write (const Matrix<T>& M, const std::string& uri) {

			if (m_iof) {
				switch (m_ios) {
					case CODRAW:  break;
					case HDF5:    return ((HDF5File*)m_iof)->Write(M,uri);
					case MATLAB:  return ((MLFile*)m_iof)->Write(M,uri);
					case ISMRM:   return ((IRDFile*)m_iof)->Write(M,uri);
					case SYNGOMR: return ((SyngoFile*)m_iof)->Write(M,uri);
					case GE:      break;
					case PHILIPS: break;
					default:      break;
				}
			}

		}

		/**
		 * @brief   Return concrete handle's status
		 */
		template<class T> Matrix<T>
		Read (const TiXmlElement* txe) const {

			if (m_iof) {
				switch (m_ios) {
					case CODRAW:  break;
					case HDF5:    return ((HDF5File*)m_iof)->Read<T>(txe);
					case MATLAB:  return ((MLFile*)m_iof)->Read<T>(txe);
					case ISMRM:   return ((IRDFile*)m_iof)->Read<T>(txe);
					case SYNGOMR: return ((SyngoFile*)m_iof)->Read<T>(txe);
					case GE:      break;
					case PHILIPS: break;
					default:      break;
				}
			}

		}

		/**
		 * @brief   Return concrete handle's status
		 */
		template <class T>
		bool Write (const Matrix<T>& M, const TiXmlElement* txe) {

			if (m_iof) {
				switch (m_ios) {
					case CODRAW:  break;
					case HDF5:    return ((HDF5File*)m_iof)->Write(M,txe);
					case MATLAB:  return ((MLFile*)m_iof)->Write(M,txe);
					case ISMRM:   return ((IRDFile*)m_iof)->Write(M,txe);
					case SYNGOMR: return ((SyngoFile*)m_iof)->Write(M,txe);
					case GE:      break;
					case PHILIPS: break;
					default:      break;
				}
			}

		}



		/**
		 * @brief   Return concrete handle's status
		 */
		error_code CleanUp () const {

			if (m_iof)
				return m_iof->CleanUp();

			return NULL_FILE_HANDLE;

		}


		/**
		 * @brief   Return concrete handle's status
		 */
		error_code Close () const {

			if (m_iof)
				return m_iof->CleanUp();

			return NULL_FILE_HANDLE;

		}

	private:


		void Concretize (const std::string& fname, const IOMode mode,
				const Params& params, const bool verbosity) {

			switch (m_ios) {
				case CODRAW:  break;
				case HDF5:    m_iof = (IOFile*) new HDF5File  (fname, mode, params, verbosity); break;
				case MATLAB:  m_iof = (IOFile*) new MLFile    (fname, mode, params, verbosity); break;
				case ISMRM:   m_iof = (IOFile*) new IRDFile   (fname, mode, params, verbosity); break;
				case SYNGOMR: m_iof = (IOFile*) new SyngoFile (fname, mode, params, verbosity); break;
				case GE:      break;
				case PHILIPS: break;
				default:      printf ("Failed to make IO context concrete!\n ");  break;
			}

		}

		IOStrategy m_ios;
		IOFile* m_iof; /**< @brief  My actual context */

	};

	/**
	 * @brief         Convenience interface to IOContext class. Close to MATLAB behaviour of fopen.
	 *
	 * @param  fname  File name (full path)
	 * @param  mode   Access mode ("r" || "w")
	 *
	 * @return        File IO context for further use.
	 */
	inline static
	IOContext fopen (const std::string& fname, const std::string& mode) {
		return IOContext (fname, mode);
	}

	/**
	 * @brief         Convenience interface to IOContext class. Close to MATLAB behaviour of fwrite.
	 *
	 * @param  f      IO context, which has been created with fopen.
	 * @param  M      Data matrix for writing.
	 *
	 * @return        Number of written elements
	 */
	template<class T> inline static
	size_t fwrite (IOContext& f, const Matrix<T>& M) {
		f.Write (M, M.GetClassName());
	}

	/**
	 * @brief         Convenience interface to IOContext class. Close to MATLAB behaviour of fread.
	 *
	 * @param  f      IO context, which has been created with fopen.
	 * @param  name   Name of data object, which needs to be read from the file.
	 *
	 * @return        Data matrix
	 */
	template<class T> inline static
	Matrix<T> fread (const IOContext& f, const std::string& name) {
		return f.Read<T> (name);
	}

	/**
	 * @brief         Convenience interface to IOContext class. Close to MATLAB behaviour of fclose.
	 *
	 * @param  f      IO context, which has been created with fopen.
	 *
	 * @return        Success (Any value > 0 indicates a problem).
	 */
	static int
	fclose (IOContext& f) {
		return (int)f.Close();
	}



}}}


#endif /* __IOCONTEXT_HPP__ */
