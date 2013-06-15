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

#include "HDF5File.hpp"
#include "SyngoFile.hpp"
#include "CODFile.hpp"
#include "Demangle.hpp"

#ifdef HAVE_ISMRMRD_HDF5_H
#include "ISMRMRD.hpp"
#endif

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
	enum IOStrategy {CODRAW = 0, HDF5, MATLAB, ISMRM, NIFTI, SYNGO, GE, PHILIPS, NO_STRATEGY};

	template<IOStrategy T> struct IOTraits;
	
	template<> struct IOTraits<HDF5> {
		typedef HDF5File IOClass;
		static const
		std::string Suffix () {
			return ".h5";
		}
		static const
		std::string CName () {
			return "hdf5";
		}
		inline static IOFile*
        Open (const std::string& fname, const IOMode mode,
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
    
	template<>
	struct IOTraits<ISMRM> {
		typedef HDF5File IOClass;

		static const std::string 
		Suffix () {
			return ".ird";
		}
		static const std::string 
		CName () {
			return "ismrm";
		}

#ifdef HAVE_ISMRMRD_HDF5_H        
		inline static IOFile*
        Open (const std::string& fname, const IOMode mode,
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
#endif        
	};
	
	template<>
	struct IOTraits<MATLAB> {

		static const
		std::string Suffix () {
			return ".mat";
		}
		static const
		std::string CName () {
			return "matlab";
		}

#ifdef HAVE_MAT_H
		typedef MLFile IOClass;
        
		inline static IOFile*
        Open (const std::string& fname, const IOMode mode,
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
#endif
	};
    
	template<>
	struct IOTraits<CODRAW> {
		typedef CODFile IOClass;
        
		static const
		std::string Suffix () {
			return ".cod";
		}
		static const
		std::string CName () {
			return "codraw";
		}
		inline static IOFile*
        Open (const std::string& fname, const IOMode mode,
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
    
	template<>
	struct IOTraits<NIFTI> {
		
		static const std::string 
		Suffix () {
			return ".nii";
		}
		static const std::string 
		CName () {
			return "nifti";
		}

#ifdef HAVE_NIFTI1_IO_H
		typedef NIFile IOClass;

		inline static IOFile*
		Open (const std::string& fname, const IOMode mode,
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
#endif
	};

	
	template<>
	struct IOTraits<SYNGO> {
		typedef SyngoFile IOClass;
		
		static const std::string
		Suffix () {
			return ".dat";
		}
		static const std::string
		CName () {
			return "syngo";
		}
		inline static IOFile*
        Open (const std::string& fname, const IOMode mode,
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
		template<class T> inline static bool
		Write (IOFile*, const Matrix<T>&, const std::string&) {
			return false;
		}
		template<class T> inline static bool
		Write (IOFile*, const Matrix<T>&, const TiXmlElement*) {
			return false;
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
			IOMode mode = (fmode.compare("r") == 0 || fmode.compare("rb") == 0) ? READ : WRITE;

			this->Concretize(fname, mode, Params(), false);

		}



		/**
		 *
		 */
		IOContext (const std::string& fname, const IOStrategy& ios = HDF5,
				const IOMode& mode = READ, const Params& params = Params(),
				const bool verbosity = false) :
			m_iof(0) , m_ios(ios) {

			this->Concretize(fname, mode, params, verbosity);

		}


		/**
		 *
		 */
		IOContext (const TiXmlElement* txe, const std::string& base = ".",
				const IOMode mode = READ, const bool verbosity = false) :
			m_iof(0), m_ios(HDF5) {

			if (!txe)
				printf ("Ouch, XML element for construction is 0!\n");

			std::string fname = base + "/" + std::string(txe->Attribute("fname"));
			std::string ftype (txe->Attribute("ftype"));

			m_ios = TypeByName (ftype);

			this->Concretize(fname, mode, Params(), verbosity);

		}


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
				switch (m_ios) {
				case CODRAW: break;
				case HDF5:   return IOTraits<  HDF5>::Read<T>(m_iof, uri);
#ifdef HAVE_MAT_H
				case MATLAB: return IOTraits<MATLAB>::Read<T>(m_iof, uri);
#endif
#ifdef HAVE_ISMRMRD_HDF5_H
				case ISMRM:  return IOTraits< ISMRM>::Read<T>(m_iof, uri);
#endif
#ifdef HAVE_NIFTI1_IO_H
				case NIFTI:  return IOTraits< NIFTI>::Read<T>(m_iof, uri);
#endif
				case SYNGO:  return IOTraits< SYNGO>::Read<T>(m_iof, uri);
				case GE:     break;
				case PHILIPS: break;
				default:     break;
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
				case CODRAW: break;
				case HDF5:   return IOTraits<  HDF5>::Write<T>(m_iof, M, uri);
#ifdef HAVE_MAT_H
				case MATLAB: return IOTraits<MATLAB>::Write<T>(m_iof, M, uri);
#endif
#ifdef HAVE_ISMRMRD_HDF5_H
				case ISMRM:  return IOTraits< ISMRM>::Write<T>(m_iof, M, uri);
#endif
#ifdef HAVE_NIFTI1_IO_H
				case NIFTI:  return IOTraits< NIFTI>::Write<T>(m_iof, M, uri);
#endif
				case SYNGO:  return IOTraits< SYNGO>::Write<T>(m_iof, M, uri);
				case GE:     break;
				case PHILIPS: break;
				default:     break;
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
				case CODRAW: break;
				case HDF5:   return IOTraits<  HDF5>::Read<T>(m_iof, txe);
#ifdef HAVE_MAT_H
				case MATLAB: return IOTraits<MATLAB>::Read<T>(m_iof, txe);
#endif
#ifdef HAVE_ISMRMRD_HDF5_H
				case ISMRM:  return IOTraits< ISMRM>::Read<T>(m_iof, txe);
#endif
#ifdef HAVE_NIFTI1_IO_H
				case NIFTI:  return IOTraits< NIFTI>::Read<T>(m_iof, txe);
#endif
				case SYNGO:  return IOTraits< SYNGO>::Read<T>(m_iof, txe);
				case GE:     break;
				case PHILIPS: break;
				default:     break;
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
				case CODRAW: break;
				case HDF5:   return IOTraits<  HDF5>::Write<T>(m_iof, M, txe);
#ifdef HAVE_MAT_H
				case MATLAB: return IOTraits<MATLAB>::Write<T>(m_iof, M, txe);
#endif
#ifdef HAVE_ISMRMRD_HDF5_H
				case ISMRM:  return IOTraits< ISMRM>::Write<T>(m_iof, M, txe);
#endif
#ifdef HAVE_NIFTI1_IO_H
				case NIFTI:  return IOTraits< NIFTI>::Write<T>(m_iof, M, txe);
#endif
				case SYNGO:  return IOTraits< SYNGO>::Write<T>(m_iof, M, txe);
				case GE:     break;
				case PHILIPS: break;
				default:     break;
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
            
			std::string lfname = strtolower (fname);

			if      (HasSuffix (lfname, IOTraits<HDF5>::Suffix()))
				return HDF5;
			else if (HasSuffix (lfname, IOTraits<MATLAB>::Suffix()))
				return MATLAB;
			else if (HasSuffix (lfname, IOTraits<CODRAW>::Suffix()))
				return CODRAW;
			else if (HasSuffix (lfname, IOTraits<NIFTI>::Suffix()))
				return NIFTI;
			else if (HasSuffix (lfname, IOTraits<SYNGO>::Suffix()))
				return SYNGO;
			else
				return HDF5;
            
		}


		IOStrategy TypeByName (const std::string& name) const {

			std::string lname = strtolower (name);

			if      (lname.compare(IOTraits<CODRAW>::CName()) == 0)
				return CODRAW;
			else if (lname.compare(IOTraits<MATLAB>::CName()) == 0)
				return MATLAB;
			else if (lname.compare(IOTraits<HDF5>::CName())   == 0)
				return HDF5;
			else if (lname.compare(IOTraits<MATLAB>::CName()) == 0)
				return NIFTI;
			else if (lname.compare(IOTraits<MATLAB>::CName()) == 0)
				return SYNGO;
			else
				return HDF5;

		}


		void Concretize (const std::string& fname, const IOMode mode,
			const Params& params, const bool verbosity) {

			switch (m_ios) {
				case CODRAW:  break;
				case HDF5:   m_iof = IOTraits<  HDF5>::Open(fname, mode, params, verbosity); break;
#ifdef HAVE_MAT_H
				case MATLAB: m_iof = IOTraits<MATLAB>::Open(fname, mode, params, verbosity); break;
#endif
#ifdef HAVE_ISMRMRD_HDF5_H
				case ISMRM:  m_iof = IOTraits< ISMRM>::Open(fname, mode, params, verbosity); break;
#endif
#ifdef HAVE_NIFTI1_IO_H
				case NIFTI:  m_iof = IOTraits< NIFTI>::Open(fname, mode, params, verbosity); break;
#endif
				case SYNGO:  m_iof = IOTraits< SYNGO>::Open(fname, mode, params, verbosity); break;
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
	IOContext fopen (const std::string& fname, const std::string& mode = "rb") {
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
