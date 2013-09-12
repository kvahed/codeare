#ifndef __TYPE_TRAITS_HPP__
#define __TYPE_TRAITS_HPP__

#include <complex>
#include <typeinfo>
#include <iostream>
#ifndef MSVC60
#  include <boost/format.hpp>
#endif

typedef std::complex<float>  cxfl;
typedef std::complex<double> cxdb;

template <class, class>
struct SameType { static const bool val = false; };

template <class T>
struct SameType<T,T> { static const bool val = true; };


template<class T> struct TypeTraits {
	inline static const std::type_info& Info () {
		return typeid(T);
	}
};

template<> struct TypeTraits<float> {
	typedef float  T;
	typedef float RT;
	typedef cxfl  CT;
	static const bool Supported;
	inline static const std::string Name () {
		return std::string("single");
	}
	inline static const bool IsComplex() {
		return SameType<T,CT>::val;
	}
	inline static const bool IsReal() {
		return SameType<T,RT>::val;
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
#ifndef MSVC60
        os << boost::format("%+.4e") % t;
#else
        os << t;
#endif
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
	inline static const bool IsComplex() {
		return SameType<T,CT>::val;
	}
	inline static const bool IsReal() {
		return SameType<T,RT>::val;
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
#ifndef MSVC60
        os << boost::format("%+.4e") % t;
#else
        os << t;
#endif
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
	inline static const bool IsComplex() {
		return SameType<T,CT>::val;
	}
	inline static const bool IsReal() {
		return SameType<T,RT>::val;
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
#ifndef MSVC60
        os << boost::format("%+.4e+1i*%+.4e") % t.real() % t.imag();
#else
        os << t;
#endif
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
	inline static const bool IsComplex() {
		return SameType<T,CT>::val;
	}
	inline static const bool IsReal() {
		return SameType<T,RT>::val;
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
#ifndef MSVC60
        os << boost::format("%+.4e+1i*%+.4e") % t.real() % t.imag();
#else
        os << t;
#endif
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
	inline static const bool IsComplex() {
		return SameType<T,CT>::val;
	}
	inline static const bool IsReal() {
		return SameType<T,RT>::val;
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
#ifndef MSVC60
        os << boost::format("%+.4e") % t;
#else
        os << t;
#endif
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
	inline static const bool IsComplex() {
		return SameType<T,CT>::val;
	}
	inline static const bool IsReal() {
		return SameType<T,RT>::val;
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
#ifndef MSVC60
        os << boost::format("%d") % (unsigned short)t;
#else
        os << (unsigned short)t;
#endif
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
	inline static const bool IsComplex() {
		return SameType<T,CT>::val;
	}
	inline static const bool IsReal() {
		return SameType<T,RT>::val;
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
#ifndef MSVC60
        os << boost::format("%li") % t;
#else
        os << t;
#endif
        return os;
    }
};

template<> struct TypeTraits<int> {
	typedef long   T;
	typedef long   RT;
	typedef void   CT;
	static const bool Supported;
	inline static const std::string Name () {
		return std::string("int");
	}
	inline static const bool IsComplex() {
		return SameType<T,CT>::val;
	}
	inline static const bool IsReal() {
		return SameType<T,RT>::val;
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
#ifndef MSVC60
        os << boost::format("%d") % t;
#else
        os << t;
#endif
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
	inline static const bool IsComplex() {
		return SameType<T,CT>::val;
	}
	inline static const bool IsReal() {
		return SameType<T,RT>::val;
	}
    inline static std::ostream& print (std::ostream& os, const T t) {
#ifndef MSVC60
        os << (unsigned long)t;
#else
        os << t;
#endif
        return os;
    }
};


#endif // __TYPE_TRAITS_HPP__
