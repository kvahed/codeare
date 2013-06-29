/*
 * Params.hpp
 *
 *  Created on: Mar 1, 2013
 *      Author: kvahed
 */

#ifndef PARAMS_HPP_
#define PARAMS_HPP_

#include "Toolbox.hpp"
#include "Demangle.hpp"

#include <boost/any.hpp>
#include <assert.h>
#include <map>
#include <string>

/**
 * @brief  General parameter container
 */
class Params {



public:


	/**
	 * @brief Default constructor
	 */
	Params() {};


	/**
	 * @brief Default destructor
	 */
	~Params() {};


	/**
	 * @brief  Parameter entry
	 */
	typedef typename std::pair<std::string, boost::any> param;


	/**
	 * @brief  Parameter container
	 */
	typedef typename std::map<std::string, boost::any> plist;


	/**
	 * @brief  Access entry (lhs)
	 */
	inline boost::any& operator [](const std::string& key) {
		plist::iterator pli = pl.find(key);
		if (pl.find(key) == pl.end()) {
			pl.insert(param(key, Toolbox::Instance()->void_any));
            pli = pl.find(key);
        }
		return pli->second;
	}


	/**
	 * @brief  Access entry (rhs)
	 */
	inline boost::any operator [](const std::string& key) const {
		plist::const_iterator pi = pl.find(key);
		if (pi != pl.end())
			return pi->second;
		else
			printf ("**WARNING**: operator[] failed for key %s\n", key.c_str());
		return key;
	}



	/**
	 * @brief       Get casted parameter value
	 *
	 * @param  key  Key
	 * @return      Casted value
	 */
	template <class T> inline
	T Get (const std::string& key) const {
		const boost::any& ba = (*this)[key];
		try {
			return boost::any_cast<T>(ba);
		} catch (const boost::bad_any_cast& e) {
			printf ("**WARNING**: Failed to retrieve %s - %s.\n             Requested %s - have %s.\n",
					key.c_str(), e.what(),
                    demangle(typeid(T).name()).c_str(),
                    demangle(ba.type().name()).c_str());
		}
		return T(0);
	}


	/**
	 * @brief       Get casted parameter value
	 *
	 * @param  key  Key
	 * @return      Casted value
	 */
	template <class T> inline T
	Get (const char* key) const {
		return Get<T> (std::string(key));
	}


	/**
	 * @brief       Get casted parameter value
	 *
	 * @param  key  Key
	 * @param  val  Casted value
	 */
	template <class T>
	inline void
	Get (const char* key, T& val) const {
		val = Get<T>(key);
	}

	/**
	 * @brief       Do we have a parameter for key
	 *
	 * @param  key  Key
	 * @return      True if yes.
	 */
	bool exists (const std::string& key) const {
		return (pl.find(key) != pl.end());
	}


	/**
	 * @brief       Overwrite/Add value to parameters
	 *
	 * @param  key  Key
	 * @param  val  Value
	 */
	inline void Set (const std::string& key, const boost::any& val) {
		if (pl.find(key) != pl.end())
			pl.erase(key);
		pl.insert(param(key, val));
	}


	/**
	 * @brief       Overwrite/Add value to parameters
	 *
	 * @param  key  Key
	 * @param  val  Value
	 */
	inline void Set (const char* key, const boost::any& val) {
		Set (std::string(key), val);
	}


	/**
	 * @brief        Get string representation of mapping
	 *
	 * @return       String representation of workspace content
	 */
	const char*
	c_str            () const {

		std::string sb;

		for (plist::const_iterator i = pl.begin(); i != pl.end(); i++)
			sb += i->first + "\t" + demangle(i->second.type().name()).c_str() + "\n";

		return sb.c_str();

	}



private:


	plist pl; /**< @brief Parameter list */

};


/**
 * @brief            Dump to ostream
 *
 * @param  os        Output stream
 * @param  w         Workspace
 * @return           The output stream
 */
inline static std::ostream&
operator<< (std::ostream& os, const Params& p) {
	os << p.c_str();
	return os;
}

#endif /* PARAMS_HPP_ */
