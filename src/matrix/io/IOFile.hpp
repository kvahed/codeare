#ifndef __IOFILE_HPP__
#define __IOFILE_HPP__

#include "Matrix.hpp"
#include "Complex.hpp"
#include "Params.hpp"
#include "Algos.hpp"
#include "Print.hpp"

#include "tinyxml/tinyxml.h"
#include "tinyxml/xpath_static.h"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace codeare {
namespace matrix  {
namespace io      {

	enum dtype {RLFL, RLDB, CXFL, CXDB, LONG, SHRT};
	enum IOMode {READ, WRITE};

	static const char* IOModeStr[] {"READ_ONLY", "READ_WRITE"};

	enum FileIOException {FILE_NOT_FOUND, OPEN_RO_FAILED, OPEN_RW_FAILED, UNKNOWN_SYNGO_DATATYPE};

    
	/**
	 * @brief           Verbose wrapper around fwrite
	 *
	 * @param  d        Data repository to write from
	 * @param  sz       Size of individual elements
	 * @param  n        Number of elements
	 * @param  f        File handle
	 * @param  desc     Description
	 * @return          Success
	 */
	template<class T> inline static bool mwrite (const T* d, const size_t n, FILE* f, 
		std::string desc) {
		size_t sz = sizeof(T);
		if (size_t l = fwrite (d, sz, n, f) != n) {
			printf("File write error - %s: %zu != %zu!\n", desc.c_str(), l, n);
			return false;
		}
		return true;
	}

	/**
	 * @brief           Verbose wrapper around fwrite
	 *
	 * @param  d        Data repository to write from
	 * @param  sz       Size of individual elements
	 * @param  n        Number of elements
	 * @param  f        File handle
	 * @param  desc     Description
	 * @return          Success
	 */
	template<class T> inline static bool
	mwrite (const std::vector<T>& d, FILE* f, std::string desc) {
		size_t n = d.size();
		size_t sz = sizeof(T);
		if (size_t l = fwrite (&d.front(), sz, n, f) != n) {
			printf("File write error - %s: %li != %li!\n", desc.c_str(), l, n);
			return false;
		}
		return true;
	}

	/**
	 * @brief           Verbose wrapper around fwrite
	 *
	 * @param  d        Data repository to write from
	 * @param  sz       Size of individual elements
	 * @param  n        Number of elements
	 * @param  f        File handle
	 * @param  desc     Description
	 * @return          Success
	 */
	template<class T> inline static bool
	mwrite (const Vector<T>& d, FILE* f, std::string desc) {
		size_t n = d.size();
		size_t sz = sizeof(T);
		if (size_t l = fwrite (d.ptr(), sz, n, f) != n) {
			printf("File write error - %s: %zu != %zu!\n", desc.c_str(), l, n);
			return false;
		}
		return true;
	}



	/**
	 * @brief           Verbose wrapper around fread
	 *
	 * @param  d        Place holder to read into
	 * @param  sz       Size of individual elements
	 * @param  n        Number of elements
	 * @param  f        File handle
	 * @param  desc     Description
	 * @return          Success
	 */
	template<class T> inline static bool mread (T* d, const size_t n, FILE* f, 
		const std::string desc) {
		size_t sz = sizeof(T);
		if (size_t l = fread (d, sz, n, f) != n) {
			printf("File read error - %s: %zu != %zu!\n", desc.c_str(), l, n);
			return false;
		}
		return true;
	}


	/**
	 * @brief           Verbose wrapper around fread
	 *
	 * @param  d        Place holder to read into
	 * @param  sz       Size of individual elements
	 * @param  n        Number of elements
	 * @param  f        File handle
	 * @param  desc     Description
	 * @return          Success
	 */
	template<class T> inline static bool
	mread (std::vector<T>& d, FILE* f, const std::string desc) {

		size_t sz = sizeof(T);
		size_t n  = d.size();

		if (size_t l = fread (&d[0], sz, n, f) != n) {
			printf("File read error - %s: %li != %li!\n", desc.c_str(), l, n);
			return false;
		}

		return true;

	}

	/**
	 * @brief           Verbose wrapper around fread
	 *
	 * @param  d        Place holder to read into
	 * @param  sz       Size of individual elements
	 * @param  n        Number of elements
	 * @param  f        File handle
	 * @param  desc     Description
	 * @return          Success
	 */
	template<class T> inline static bool
	mread (Vector<T>& d, FILE* f, const std::string desc) {

		size_t sz = sizeof(T);
		size_t n  = d.size();

		if (size_t l = fread (&d[0], sz, n, f) != n) {
			printf("File read error - %s: %zu != %zu!\n", desc.c_str(), l, n);
			return false;
		}

		return true;

	}


    /**
     * brief          Does a file exist?
     * 
     * @param  fname  File name
     * @return        Does it exist?
     */
    inline static bool
    fexists (const char* fname) {
    	assert (fname[0] != '\0');
        ifstream fs (fname);
		return (fs.is_open());
    }
    

    /**
     * brief          Does a file exist?
     *
     * @param  fname  File name
     * @return        Does it exist?
     */
    inline static bool
    fexists (const std::string& fname) {
        return (fexists(fname.c_str()));
    }




	// TODO: deliver some crap
	template<class T>
	inline static bool is_mat (boost::any oper) {

		try {
			boost::any_cast<boost::shared_ptr<Matrix<T> > >(oper);
		} catch(const boost::bad_any_cast &) {
			return false;
		}

		return true;

	}

	/**
	 * @brief     Base class for file IO
	 */
	class IOFile {

	public:

		typedef std::map<std::string, boost::any> DataStore;
		typedef std::pair<std::string, boost::any> DataEntry;


		/**
		 * @brief Construct with file
		 *
		 * @param  fname    Filename
		 * @param  mode     IO mode (READ/WRITE)
		 * @param  params   Optional parameters
		 * @param  verbose  Verbose output? (default: true)
		 */
		IOFile (const std::string& fname, const IOMode& mode = READ,
				Params params = Params(), const bool& verbose = true) :
			m_fname(fname),
			m_status(OK),
			m_params(params),
			m_fopen(false),
			m_verb(verbose),
			m_alloc(false),
			m_initialised (false),
			m_mode (mode) {

			// TODO: Check if file exists, set status

		}

		/**
		 * @brief  Default destructor
		 */
		virtual ~IOFile () {
			CleanUp();
		}


		/**
		 * @brief Clean up
		 */
		virtual error_code
		CleanUp () {
			return OK;
		}


		/**
		 * @brief  File name
		 *
		 * @return  File name
		 */
		inline std::string
		FileName() const {
			return m_fname;
		}


		/**
		 * @brief  Verbosity
		 *
		 * @return Verbosity
		 */
		inline bool
		Verbosity() const {
			return m_verb;
		}


		/**
		 * @brief  Status
		 *
		 * @return Status
		 */
		inline error_code
		Status () const {
			return m_status;
		}


		/**
		 * @brief  Memory allocated?
		 *
		 * @return Memory allocated?
		 */
		inline bool
		Allocated () const {
			return m_alloc;
		}


		/**
		 * @brief  Holding write lock?
		 *
		 * @return Write lock?
		 */
		inline bool
		Locked () const {
			return (m_mode == WRITE);
		}


		/**
		 * @brief Read a particular data set from file
		 *
		 * @return  Success
		 */
		template<class T> Matrix<T>
		Read (const std::string& uri) const;


		/**
		 * @brief  Write data to file
		 *
		 * @return  Success
		 */
		template<class T> bool
		Write (const Matrix<T>& M, const std::string& uri);


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

	protected:

		std::string m_fname; /**< @brief  File name */

		error_code  m_status; /**< @brief Status */

		IOMode      m_mode;   /**< @brief IO mode */

		Params      m_params; /**< @brief Additional parameters */
		bool        m_fopen;  /**< @brief File open? */
		bool        m_verb;   /**< @brief Verbosity */
		bool        m_alloc;  /**< @brief Memory allocated? */
		bool        m_initialised; /**< @brief Initialised */


		DataStore   m_data; /**< Data */


		/**
		 * @brief Default constructor
		 */
		IOFile () :
			m_fname(""),
			m_status(OK),
			m_fopen(false),
			m_verb(true),
			m_alloc(false),
			m_initialised (false),
			m_mode (READ),
			m_params(Params()){}

		/**
		 * @brief Copy constructor
		 *
		 * @param  iob  To copy
		 */
		IOFile (const IOFile& iob):
			m_fname(iob.FileName()),
			m_status(iob.Status()),
			m_verb (iob.Verbosity()),
			m_alloc (iob.Allocated()),
			m_fopen (false),
			m_initialised (false),
			m_mode (READ),
			m_params (Params()) {}

	};

} // namespace io
} // namespace matrix
} // namespace codeare

#endif
