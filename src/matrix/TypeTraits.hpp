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
    inline static const std::string Abbrev () {
        return std::string("float");
    }
    inline static bool IsComplex() {
        return SameType<T,CT>::val;
    }
    inline static bool IsReal() {
        return SameType<T,RT>::val;
    }
    inline static RT Real (const T t) {return t;}
    inline static RT Imag (const T) {return 0.;}
    inline static RT Abs (const T t) {return t;}
    inline static RT Arg (const T) {return 0.;}
    inline static RT Conj (const T t) {return t;}
    inline static RT Pow (const T& x, const T& y) {return std::pow(x,y);}
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
    inline static const std::string Abbrev () {
        return std::string("double");
    }
    inline static bool IsComplex() {
        return SameType<T,CT>::val;
    }
    inline static bool IsReal() {
        return SameType<T,RT>::val;
    }
    inline static RT Real (const T t) {return t;}
    inline static RT Imag (const T) {return 0.;}
    inline static RT Abs (const T t) {return t;}
    inline static RT Arg (const T) {return 0.;}
    inline static RT Conj (const T t) {return t;}
    inline static RT Pow (const T& x, const T& y) {return std::pow(x,y);}
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
    inline static const std::string Abbrev () {
        return std::string("cxfl");
    }
    inline static bool IsComplex() {
        return SameType<T,CT>::val;
    }
    inline static bool IsReal() {
        return SameType<T,RT>::val;
    }
    inline static RT Real (const T t) {return std::real(t);}
    inline static RT Imag (const T t) {return std::imag(t);}
    inline static RT Abs (const T t) {return std::abs(t);}
    inline static RT Arg (const T t) {return std::arg(t);}
    inline static CT Conj (const T t) {return std::conj(t);}
    inline static CT Pow (const CT& x, const T& y) {return pow(x,y);}
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
    inline static const std::string Abbrev () {
        return std::string("cxdb");
    }
    inline static bool IsComplex() {
        return SameType<T,CT>::val;
    }
    inline static bool IsReal() {
        return SameType<T,RT>::val;
    }
    inline static RT Real (const T t) {return std::real(t);}
    inline static RT Imag (const T t) {return std::imag(t);}
    inline static RT Abs (const T t) {return std::abs(t);}
    inline static RT Arg (const T t) {return std::arg(t);}
    inline static CT Conj (const T t) {return std::conj(t);}
    inline static CT Pow (const CT& x, const T& y) {return pow(x,y);}
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
    inline static const std::string Abbrev () {
        return std::string("short");
    }
    inline static bool IsComplex() {
        return SameType<T,CT>::val;
    }
    inline static bool IsReal() {
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
    inline static const std::string Abbrev () {
        return std::string("uchar");
    }
    inline static bool IsComplex() {
        return SameType<T,CT>::val;
    }
    inline static  bool IsReal() {
        return SameType<T,RT>::val;
    }
    inline static T Pow (const T& x, const T& y) { return x;}
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
    inline static const std::string Abbrev () {
        return std::string("long");
    }
    inline static bool IsComplex() {
        return SameType<T,CT>::val;
    }
    inline static bool IsReal() {
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
    inline static const std::string Abbrev () {
        return std::string("int");
    }
    inline static bool IsComplex() {
        return SameType<T,CT>::val;
    }
    inline static bool IsReal() {
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
    inline static const std::string Abbrev () {
        return std::string("size_t");
    }
    inline static bool IsComplex() {
        return SameType<T,CT>::val;
    }
    inline static bool IsReal() {
        return SameType<T,RT>::val;
    }
    inline static T Pow (const T& x, const T& y) { return (T)pow((double)x,(double)y);} 
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
