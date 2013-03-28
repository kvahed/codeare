#ifndef __IOFILE_HPP__
#define __IOFILE_HPP__

#include "Matrix.hpp"
#include "Params.hpp"

namespace codeare {
namespace matrix  {
namespace io      {

	enum IOMode {READ, WRITE};

	// TODO: deliver some crap
	template<class T>
	inline static bool is_mat (boost::any oper) {

		try {
			boost::any_cast<Ptr<Matrix<T> > >(oper);
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

			while (!m_data.empty()) {
				delete boost::any_cast<Ptr<Matrix<cxfl> > >(m_data.begin()->second);
				m_data.erase(m_data.begin());
			}

			return OK;

		}


		/**
		 * @brief  File name
		 *
		 * @return  File name
		 */
		std::string
		FileName() const {
			return m_fname;
		}


		/**
		 * @brief  Verbosity
		 *
		 * @return Verbosity
		 */
		bool
		Verbosity() const {
			return m_verb;
		}


		/**
		 * @brief  Status
		 *
		 * @return Status
		 */
		error_code
		Status () const {
			return m_status;
		}


		/**
		 * @brief  Memory allocated?
		 *
		 * @return Memory allocated?
		 */
		bool
		Allocated () const {
			return m_alloc;
		}


		/**
		 * @brief  Holding write lock?
		 *
		 * @return Write lock?
		 */
		bool
		Locked () const {
			return (m_mode == WRITE);
		}


		/**
		 * @brief Read a particular data set from file
		 *
		 * @return  Success
		 */
		virtual error_code
		Read (const std::string& dname = "") = 0;


		/**
		 * @brief  Write data to file
		 *
		 * @return  Success
		 */
		virtual error_code
		Write (const std::string& dname = "") = 0;

		/**
		 * @brief  Get data by name
		 *
		 * @return  Any data type
		 */
		boost::any&
		GetAny (const std::string& dname) {

			std::map<std::string, boost::any>::iterator it;
			it = m_data.find(dname);

			if (it == m_data.end()) {
				printf ("  ERROR: No dataset by the name %s found?!\n",
						dname.c_str());
				return Toolbox::Instance()->void_any;
			}

			return it->second;

		}


		/**
		 * @brief  Get casted Matrix by name
		 *
		 * @return  Any data type
		 */
		template<class T> Matrix<T>&
		Get (const std::string& dname) {
			return *boost::any_cast<Ptr<Matrix<T> > >(GetAny(dname));
		}





	protected:

		std::string m_fname; /**< @brief  File name */

		error_code  m_status; /**< @brief Status */

		IOMode      m_mode;   /**< @Brief IO mode */

		bool        m_fopen;  /**< @brief File open? */
		bool        m_verb;   /**< @brief Verbosity */
		bool        m_alloc;  /**< @brief Memory allocated? */
		bool        m_initialised; /**< @brief Initialised */

		Params      m_params; /**< @brief Additional parameters */

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
