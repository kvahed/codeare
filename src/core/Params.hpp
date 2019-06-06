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

inline static int signed_cast (const boost::any& b)  {
    int ret;
    try {
        ret = boost::any_cast<int> (b);
    } catch (const boost::bad_any_cast&) {
        try {
            ret = boost::any_cast<long> (b);
        } catch (const boost::bad_any_cast&) {
            try {
                ret = boost::any_cast<short> (b);
            } catch (const boost::bad_any_cast&) {
                try {
                    ret = boost::any_cast<size_t> (b);
                } catch (const boost::bad_any_cast&) {
                    try {
                        ret = boost::any_cast<unsigned> (b);
                    } catch (const boost::bad_any_cast&) {
                        try {
                            ret = boost::any_cast<unsigned long> (b);
                        } catch (const boost::bad_any_cast& e) {
                            throw e;
                        }
                    }
                }
            }
        }
    }
    return ret;
}

inline static size_t unsigned_cast (const boost::any& b) {
    size_t ret;
    try {
        ret = boost::any_cast<int> (b);
    } catch (const boost::bad_any_cast&) {
        try {
            ret = boost::any_cast<long> (b);
        } catch (const boost::bad_any_cast&) {
            try {
                ret = boost::any_cast<short> (b);
            } catch (const boost::bad_any_cast&) {
                try {
                    ret = boost::any_cast<size_t> (b);
                } catch (const boost::bad_any_cast&) {
                    try {
                        ret = boost::any_cast<unsigned> (b);
                    } catch (const boost::bad_any_cast&) {
                        try {
                            ret = boost::any_cast<unsigned long> (b);
                        } catch (const boost::bad_any_cast& e) {
                            throw e;
                        }
                    }
                }
            }
        }
    }
    return ret;
}

inline static std::complex<float> complex_cast (const boost::any& b)  {
    std::complex<float> ret;
    try {
        ret = boost::any_cast<std::complex<float> > (b);
    } catch (const boost::bad_any_cast&) {
        try {
            ret = boost::any_cast<std::complex<double> > (b);
        } catch (const boost::bad_any_cast&) {
            try {
                ret = boost::any_cast<float> (b);
            } catch (const boost::bad_any_cast&) {
                try {
                    ret = boost::any_cast<double> (b);
                } catch (const boost::bad_any_cast&) {
                    try {
                        ret = boost::any_cast<int> (b);
                    } catch (const boost::bad_any_cast&) {
                        try {
                            ret = boost::any_cast<long> (b);
                        } catch (const boost::bad_any_cast&) {
                            try {
                                ret = boost::any_cast<short> (b);
                            } catch (const boost::bad_any_cast&) {
                                try {
                                    ret = boost::any_cast<size_t> (b);
                                } catch (const boost::bad_any_cast&) {
                                    try {
                                        ret = boost::any_cast<unsigned> (b);
                                    } catch (const boost::bad_any_cast&) {
                                        try {
                                            ret = boost::any_cast<unsigned long> (b);
                                        } catch (const boost::bad_any_cast& e) {
                                            throw e;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return ret;
}

inline static float fp_cast (const boost::any& b) {
    float ret;
    try {
        ret = boost::any_cast<float> (b);
    } catch (const boost::bad_any_cast&) {
        try {
            ret = boost::any_cast<double> (b);
        } catch (const boost::bad_any_cast&) {
            try {
                ret = boost::any_cast<int> (b);
            } catch (const boost::bad_any_cast&) {
                try {
                    ret = boost::any_cast<long> (b);
                } catch (const boost::bad_any_cast&) {
                    try {
                        ret = boost::any_cast<short> (b);
                    } catch (const boost::bad_any_cast&) {
                        try {
                            ret = boost::any_cast<size_t> (b);
                        } catch (const boost::bad_any_cast&) {
                            try {
                                ret = boost::any_cast<unsigned> (b);
                            } catch (const boost::bad_any_cast&) {
                                try {
                                    ret = boost::any_cast<unsigned long> (b);
                                } catch (const boost::bad_any_cast& e) {
                                    throw e;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return ret;
}

enum PARAMETER_MAP_EXCEPTION {KEY_NOT_FOUND};

/**
 * @brief  General parameter container
 */
class DLLEXPORT Params {


    struct UnspecParam {};

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
	 * @brief  Access entry (lhs)
	 */
	inline boost::any& operator [](const std::string& key) {
		std::map<std::string, boost::any>::iterator pli = pl.find(key);
		if (pl.find(key) == pl.end()) {
			pl.insert(std::pair<std::string, boost::any>(key, UnspecParam()));
            pli = pl.find(key);
        }
		return pli->second;
	}


	/**
	 * @brief  Access entry (rhs)
	 */
	inline boost::any operator [](const std::string& key) const {
		std::map<std::string, boost::any>::const_iterator pi = pl.find(key);
		if (pi != pl.end())
			return pi->second;
		else
			throw KEY_NOT_FOUND;
        //printf ("**WARNING**: operator[] failed for key %s\n", key.c_str());
		return key;
	}



	/**
	 * @brief       Get casted parameter value
	 *
	 * @param  key  Key
	 * @return      Casted value
	 */
	template <class T> inline T Get (const std::string& key) const {
		const boost::any& ba = (*this)[key];
		try {
			return boost::any_cast<T>(ba);
		} catch (const boost::bad_any_cast& e) {
			printf ("**WARNING**: Failed to retrieve %s - %s.\n             Requested %s - have %s.\n",
					key.c_str(), e.what(),
                    demangle(typeid(T).name()).c_str(),
                    demangle(ba.type().name()).c_str());
			throw e;
		}
		T t = T();
		return t;
	}


	/**
	 * @brief       Get casted parameter value
	 *
	 * @param  key  Key
	 * @return      Casted value
	 */
	template <class T> inline T Get (const char* key) const {
		return Get<T> (std::string(key));
	}


	/**
	 * @brief       Get casted parameter value
	 *
	 * @param  key  Key
	 * @param  val  Casted value
	 */
	template <class T> inline void Get (const char* key, T& val) const {
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
		pl.insert(std::pair<std::string, boost::any>(key, val));
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
	void Print (std::ostream& os) const {
		typedef std::map<std::string, boost::any>::const_iterator it;
		for (it i = pl.begin(); i != pl.end(); i++) {
			const boost::any& b = i->second;
			std::string v_name  = demangle(i->second.type().name()),
					    k_name  = i->first;

			os << setw(24) << k_name << " | "
			   << setw(16) << v_name << " | "
			   << setw(32);
			if (b.type() == typeid(float))
				TypeTraits<float>::print(os,boost::any_cast<float>(b));
			else if (b.type() == typeid(double))
				TypeTraits<double>::print(os,boost::any_cast<double>(b));
			else if (b.type() == typeid(cxfl))
				TypeTraits<cxfl>::print(os,boost::any_cast<cxfl>(b));
			else if (b.type() == typeid(cxdb))
				TypeTraits<cxdb>::print(os,boost::any_cast<cxdb>(b));
			else if (b.type() == typeid(int))
				os << boost::any_cast<int>(b);
			else if (b.type() == typeid(short))
				os << boost::any_cast<short>(b);
			else if (b.type() == typeid(long))
				os << boost::any_cast<long>(b);
			else if (b.type() == typeid(size_t))
				os << boost::any_cast<size_t>(b);
			else if (b.type() == typeid(cbool))
				os << boost::any_cast<cbool>(b);
			else if (b.type() == typeid(char*))
				os << boost::any_cast<char*>(b);
			else if (b.type() == typeid(const char*))
				os << boost::any_cast<const char*>(b);
			else if (b.type() == typeid(std::string))
				os << boost::any_cast<std::string>(b).c_str();
			os << "\n";
		}
	}



private:
#ifdef _MSC_VER
#pragma warning (disable : 4251)
#endif
	std::map<std::string, boost::any> pl; /**< @brief Parameter list */
#ifdef _MSC_VER
#pragma warning (default : 4251)
  #endif
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
	p.Print(os);
	return os;
}

template<class T> inline T try_to_fetch (const Params& p, const std::string& key, const T& def) {
    try {
        return p.Get<T>(key);
    } catch (const PARAMETER_MAP_EXCEPTION&) {
        std::cout << "**WARNING**: key " << key.c_str() << " not in parameter map." << std::endl;
        std::cout << "             Returning default " << def << std::endl;
    } catch (const boost::bad_any_cast&) {
        std::cout << "**WARNING**: conversion failed on parameter " <<  key.c_str() << "." << std::endl;
        std::cout << "             Returning default " << def << std::endl;
    }
    return def;
}

#endif /* PARAMS_HPP_ */


