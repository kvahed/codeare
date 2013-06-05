/*
 *  codeare Copyright (C) 2010-2013 
 *                        Kaveh Vahedipour
 *                        Forschungszentrum Juelich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but 
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301  USA
 */

#ifndef __IOCONTEXT_HPP__
#define __IOCONTEXT_HPP__

#include "ISMRMRD.hpp"
#include "HDF5File.hpp"
#include "SyngoFile.hpp"
#include "CODFile.hpp"

#ifdef HAVE_MAT_H
#include "MLFile.hpp"
#endif

#ifdef HAVE_NIFTI1_IO_H
#include "NIFile.hpp"
#endif

/*
#include "NetCDF.hpp"
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

	template<IOStrategy T>
	struct IOTraits;

	template<>
	struct IOTraits<HDF5> {
		typedef HDF5File IOClass;

		static const std::string Suffix () {
			return ".h5";
		}
		static const std::string CName () {
			return "CODRAW";
		}
		inline IOFile* Open (const std::string& fname, const IOMode mode,
				const Params& params, const bool verbosity) {
			return (IOFile*) new IOClass (fname, mode, params, verbosity);
		}
		template <class T> inline static Matrix<T>
		Read (const IOFile* iof, const std::string& uri) {
			return ((IOClass*)iof)->Read<T>(uri);
		}
		template <class T> inline static Matrix<T>
		Read (const IOFile* iof, const TiXmlElement* txe) {
			return ((IOClass*)iof)->Read<T>(txe);
		}
		template <class T> inline static bool
		Write (IOFile* iof, const Matrix<T>& M, const std::string& uri) {
			return ((IOClass*)iof)->Write(M,uri);
		}
		template <class T> inline static bool
		Write (IOFile* iof, const Matrix<T>& M, const TiXmlElement* txe) {
			return ((IOClass*)iof)->Write(M,txe);
		}
	};

#ifdef HAVE_MAT_H
	template<>
	struct IOTraits<MATLAB> {
		typedef MLFile IOClass;

		static const std::string Suffix () {
			return ".mat";
		}
		static const std::string CName () {
			return "MATLAB";
		}
		inline IOFile* Open (const std::string& fname, const IOMode mode,
				const Params& params, const bool verbosity) {
			return (IOFile*) new IOClass (fname, mode, params, verbosity);
		}
		template <class T> inline static Matrix<T>
		Read (const IOFile* iof, const std::string& uri) {
			return ((IOClass*)iof)->Read<T>(uri);
		}
		template <class T> inline static Matrix<T>
		Read (const IOFile* iof, const TiXmlElement* txe) {
			return ((IOClass*)iof)->Read<T>(txe);
		}
		template <class T> inline static bool
		Write (IOFile* iof, const Matrix<T>& M, const std::string& uri) {
			return ((IOClass*)iof)->Write(M,uri);
		}
		template <class T> inline static bool
		Write (IOFile* iof, const Matrix<T>& M, const TiXmlElement* txe) {
			return ((IOClass*)iof)->Write(M,txe);
		}
	};
#endif

	template<>
	struct IOTraits<CODRAW> {
		typedef CODFile IOClass;

		static const std::string Suffix () {
			return ".cod";
		}
		static const std::string CName () {
			return "CODRAW";
		}
		inline IOFile* Open (const std::string& fname, const IOMode mode,
				const Params& params, const bool verbosity) {
			return (IOFile*) new IOClass (fname, mode, params, verbosity);
		}
		template <class T> inline static Matrix<T>
		Read (const IOFile* iof, const std::string& uri) {
			return ((IOClass*)iof)->Read<T>(uri);
		}
		template <class T> inline static Matrix<T>
		Read (const IOFile* iof, const TiXmlElement* txe) {
			return ((IOClass*)iof)->Read<T>(txe);
		}
		template <class T> inline static bool
		Write (IOFile* iof, const Matrix<T>& M, const std::string& uri) {
			return ((IOClass*)iof)->Write(M,uri);
		}
		template <class T> inline static bool
		Write (IOFile* iof, const Matrix<T>& M, const TiXmlElement* txe) {
			return ((IOClass*)iof)->Write(M,txe);
		}
	};

#ifdef HAVE_NIFTI1_IO_H
	template<>
	struct IOTraits<NIFTI> {
		typedef NIFile IOClass;

		static const std::string Suffix () {
			return ".nii";
		}
		static const std::string CName () {
			return "NIFTI";
		}
		inline IOFile* Open (const std::string& fname, const IOMode mode,
				const Params& params, const bool verbosity) {
			return (IOFile*) new IOClass (fname, mode, params, verbosity);
		}
		template <class T> inline static Matrix<T>
		Read (const IOFile* iof, const std::string& uri) {
			return ((IOClass*)iof)->Read<T>(uri);
		}
		template <class T> inline static Matrix<T>
		Read (const IOFile* iof, const TiXmlElement* txe) {
			return ((IOClass*)iof)->Read<T>(txe);
		}
		template <class T> inline static bool
		Write (IOFile* iof, const Matrix<T>& M, const std::string& uri) {
			return ((IOClass*)iof)->Write(M,uri);
		}
		template <class T> inline static bool
		Write (IOFile* iof, const Matrix<T>& M, const TiXmlElement* txe) {
			return ((IOClass*)iof)->Write(M,txe);
		}
	};
#endif

	template<>
	struct IOTraits<SYNGOMR> {
		typedef NIFile IOClass;

		static const std::string Suffix () {
			return ".dat";
		}
		static const std::string CName () {
			return "SYNGOMR";
		}
		inline IOFile* Open (const std::string& fname, const IOMode mode,
				const Params& params, const bool verbosity) {
			return (IOFile*) new IOClass (fname, mode, params, verbosity);
		}
		template <class T> inline static Matrix<T>
		Read (const IOFile* iof, const std::string& uri) {
			return ((IOClass*)iof)->Read<T>(uri);
		}
		template <class T> inline static Matrix<T>
		Read (const IOFile* iof, const TiXmlElement* txe) {
			return ((IOClass*)iof)->Read<T>(txe);
		}
	};


	/**
	 * @brief       Interface to concrete IO implementation
	 */
	class IOContext {

	public:

		/**
		 * @brief   Default constructor
		 */
		IOContext () :
			m_iof(0), m_ios(HDF5) {}

		/**
		 * @brief   Default copy constructor
		 */
		IOContext (const IOContext& ioc) :
			m_iof(0), m_ios(HDF5) {}


		/**
		 *
		 */
		IOContext (const std::string& fname, const std::string& fmode) :
			m_iof(0), m_ios(HDF5) {

			assert (fmode.compare("r") == 0 || fmode.compare("rb") == 0 || fmode.compare("w") == 0 || fmode.compare("wb") == 0);

			m_ios = TypeBySuffix (fname);

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

			m_ios = TypeByName (ftype);

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
#ifdef HAVE_MAT_H
					case MATLAB:  return ((MLFile*)m_iof)->Read<T>(uri);
#endif
					case ISMRM:   return ((IRDFile*)m_iof)->Read<T>(uri);
#ifdef HAVE_NIFTI1_IO_H
					case NIFTI:   return ((NIFile*)m_iof)->Read<T>(uri);
#endif
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
#ifdef HAVE_MAT_H
					case MATLAB:  return ((MLFile*)m_iof)->Write(M,uri);
#endif
					case ISMRM:   return ((IRDFile*)m_iof)->Write(M,uri);
#ifdef HAVE_NIFTI1_IO_H
					case NIFTI:   return ((NIFile*)m_iof)->Write(M,uri);
#endif
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
#ifdef HAVE_MAT_H
					case MATLAB:  return ((MLFile*)m_iof)->Read<T>(txe);
#endif
					case ISMRM:   return ((IRDFile*)m_iof)->Read<T>(txe);
#ifdef HAVE_NIFTI1_IO_H
					case NIFTI:   return ((NIFile*)m_iof)->Read<T>(txe);
#endif
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
#ifdef HAVE_MAT_H
					case MATLAB:  return ((MLFile*)m_iof)->Write(M,txe);
#endif
					case ISMRM:   return ((IRDFile*)m_iof)->Write(M,txe);
#ifdef HAVE_NIFTI1_IO_H
					case NIFTI:   return ((NIFile*)m_iof)->Write(M,txe);
#endif
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


		IOStrategy TypeBySuffix (const std::string& fname) const {

			if      (HasSuffix (fname, IOTraits<HDF5>::CName()))
				return HDF5;
			else if (HasSuffix (fname, IOTraits<MATLAB>::CName()))
				return MATLAB;
			else if (HasSuffix (fname, IOTraits<CODRAW>::CName()))
				return CODRAW;
			else if (HasSuffix (fname, IOTraits<NIFTI>::CName()))
				return NIFTI;
			else if (HasSuffix (fname, IOTraits<SYNGOMR>::CName()))
				return SYNGOMR;
			else
				return HDF5;

		}


		IOStrategy TypeByName (const std::string& name) const {

			if (name.compare("CODRAW") == 0)
				return CODRAW;
			else if (name.compare("MATLAB") == 0)
				return MATLAB;
			else if (name.compare("HDF5") == 0)
				return HDF5;
			else if (name.compare("NIFTI") == 0)
				return NIFTI;
			else if (name.compare("SYNGOMR") == 0)
				return SYNGOMR;
			else
				return HDF5;

		}


		void Concretize (const std::string& fname, const IOMode mode,
			const Params& params, const bool verbosity) {

			switch (m_ios) {
				case CODRAW:  break;
				case HDF5:    m_iof = (IOFile*) new HDF5File  (fname, mode, params, verbosity); break;
#ifdef HAVE_MAT_H
				case MATLAB:  m_iof = (IOFile*) new MLFile    (fname, mode, params, verbosity); break;
#endif
				case ISMRM:   m_iof = (IOFile*) new IRDFile   (fname, mode, params, verbosity); break;
#ifdef HAVE_NIFTI1_IO_H
				case NIFTI:   m_iof = (IOFile*) new NIFile    (fname, mode, params, verbosity); break;
#endif
				case SYNGOMR: m_iof = (IOFile*) new SyngoFile (fname, mode, params, verbosity); break;
				case GE:      break;
				case PHILIPS: break;
				default:      printf ("Failed to make IO context concrete!\n ");  break;
			}

		}


		const IOStrategy Strategy () const {
			return m_ios;
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
#define fwrite(X,Y) _fwrite(X,Y,#Y)
	template<class T> inline static
	size_t _fwrite (IOContext& f, const Matrix<T>& M, const std::string& name) {
		f.Write (M, name);
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
