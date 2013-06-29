#ifndef __TYPE_TRAITS_HPP__
#define __TYPE_TRAITS_HPP__

#include <complex>
#include <typeinfo>
#include <iostream>
#include <boost/format.hpp>

typedef std::complex<float>  cxfl;
typedef std::complex<double> cxdb;

template<class T> struct TypeTraits;

template<> struct TypeTraits<float> {
	typedef float  T;
	typedef float RT;
	typedef cxfl  CT;
	static const bool Supported;
	inline static const std::string Name () {
		return std::string("single");
	}
	inline static const std::type_info& Info () {
		return typeid(T);
	}
	inline static const bool IsComplex() {
		return (typeid(T) == typeid(CT));
	}
	inline static const bool IsReal() {
		return (typeid(T) == typeid(RT));
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
        os << boost::format("%+.4e") % t;
        return os;
    }
};

template<> struct TypeTraits<double> {
	typedef double  T;
	typedef double RT;
	typedef cxdb  CT;
	static const bool Supported;

	inline static const std::string Name () {
		return std::string("double");
	}
	inline static const std::type_info& Info () {
		return typeid(T);
	}
	inline static const bool IsComplex() {
		return (typeid(T) == typeid(CT));
	}
	inline static const bool IsReal() {
		return (typeid(T) == typeid(RT));
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
        os << boost::format("%+.4e") % t;
        return os;
    }

};

template<> struct TypeTraits<cxfl> {
	typedef cxfl T;
	typedef float RT;
	typedef cxfl  CT;
	static const bool Supported;
	inline static const std::string Name () {
		return std::string("complex single");
	}
	inline static const std::type_info& Info () {
		return typeid(T);
	}
	inline static const bool IsComplex() {
		return (typeid(T) == typeid(CT));
	}
	inline static const bool IsReal() {
		return (typeid(T) == typeid(RT));
	}
	inline static const bool Validate () {
		return true;
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
        os << boost::format("%+.4e+1i*%+.4e") % t.real() % t.imag();
        return os;
    }
};

template<> struct TypeTraits<cxdb> {
	typedef cxdb  T;
	typedef double RT;
	typedef cxdb  CT;
	static const bool Supported;
	inline static const std::string Name () {
		return std::string("complex double");
	}
	inline static const std::type_info& Info () {
		return typeid(T);
	}
	inline static const bool IsComplex() {
		return (typeid(T) == typeid(CT));
	}
	inline static const bool IsReal() {
		return (typeid(T) == typeid(RT));
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
        os << boost::format("%+.4e+1i*%+.4e") % t.real() % t.imag();
        return os;
    }
};

template<> struct TypeTraits<short> {
	typedef short  T;
	typedef short  RT;
	typedef void   CT;
	static const bool Supported;
	inline static const std::string Name () {
		return std::string("short int");
	}
	inline static const std::type_info& Info () {
		return typeid(T);
	}
	inline static const bool IsComplex() {
		return (typeid(T) == typeid(CT));
	}
	inline static const bool IsReal() {
		return (typeid(T) == typeid(RT));
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
        os << boost::format("%+.4e") % t;
        return os;
    }
};

template<> struct TypeTraits<unsigned char> {
	typedef unsigned char  T;
	typedef unsigned char  RT;
	typedef void   CT;
	static const bool Supported;
	inline static const std::string Name () {
		return std::string("codeare bool");
	}
	inline static const std::type_info& Info () {
		return typeid(T);
	}
	inline static const bool IsComplex() {
		return (typeid(T) == typeid(CT));
	}
	inline static const bool IsReal() {
		return (typeid(T) == typeid(RT));
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
        os << boost::format("%d") % (unsigned short)t;
        return os;
    }
};


template<> struct TypeTraits<long> {
	typedef long   T;
	typedef long   RT;
	typedef void   CT;
	static const bool Supported;
	inline static const std::string Name () {
		return std::string("long int");
	}
	inline static const std::type_info& Info () {
		return typeid(T);
	}
	inline static const bool IsComplex() {
		return (typeid(T) == typeid(CT));
	}
	inline static const bool IsReal() {
		return (typeid(T) == typeid(RT));
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
        os << boost::format("%li") % t;
        return os;
    }
};

template<> struct TypeTraits<size_t> {
	typedef size_t  T;
	typedef size_t RT;
	typedef void   CT;
	static const bool Supported;
	inline static const std::string Name () {
		return std::string("size type");
	}
	inline static const std::type_info& Info () {
		return typeid(T);
	}
	inline static const bool IsComplex() {
		return (typeid(T) == typeid(CT));
	}
	inline static const bool IsReal() {
		return (typeid(T) == typeid(RT));
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
        os << boost::format("%li") % t;
        return os;
    }
};


#endif // __TYPE_TRAITS_HPP__
