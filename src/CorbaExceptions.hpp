/**
 * @brief Corba communication exception
 */
class DS_ServerConnectionException {
public:
	/**
	 * @brief Handle communication exceptions
	 */
	DS_ServerConnectionException () { 
		std::cerr << "CORBA server connection exception" << std::endl; 
	};
};


/**
 * @brief Corba system exception
 */

class DS_SystemException           {
public:
	/**
	 * @brief Handle system exceptions
	 */
	DS_SystemException           () { 
		std::cerr << "CORBA system exception" << std::endl; 
	};
};


/**
 * @brief Corba fatal exception
 */
class DS_FatalException            {
public:
	/**
	 * @brief Handle fatal exception
	 */
	DS_FatalException            () { 
		std::cerr << "CORBA fatal exception" << std::endl; 
	};
};


/**
 * @brief Corba unspecified exception
 */
class DS_Exception                 {
public:
	/**
	 * @brief Handle unspecified exception
	 */
	DS_Exception                 () { 
		std::cerr << "CORBA unspecified exception" << std::endl; 
	};
};

