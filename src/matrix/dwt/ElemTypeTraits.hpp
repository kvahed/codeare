# ifndef __ELEM_TYPE_TRAITS_HPP__

# define __ELEM_TYPE_TRAITS_HPP__


/******************************
 ** struct: elem_type_traits **
 **     (base struct)        **
 ******************************/
template <class T>
struct elem_type_traits
{

        /* -- */

}; // struct elem_type_traits <T>


/******************************
 ** struct: elem_type_traits **
 **     (spec: float)        **
 ******************************/
template <>
struct elem_type_traits <float>
{

    public:

        typedef float elem_type;
        typedef float value_type;

        static inline
        const char *
        print_elem_type       ( )
        {
            return "float";
        }


}; // struct elem_type_traits <float>


/******************************
 ** struct: elem_type_traits **
 **     (spec: cxfl)         **
 ******************************/
template <>
struct elem_type_traits <cxfl>
{

    public:

        typedef cxfl elem_type;
        typedef float value_type;

        static inline
        const char *
        print_elem_type       ( )
        {
            return "cxfl";
        }

}; // struct elem_type_traits <cxfl>


/******************************
 ** struct: elem_type_traits **
 **     (spec: double)       **
 ******************************/
template <>
struct elem_type_traits <double>
{

    public:

        typedef double elem_type;
        typedef double value_type;

        static inline
        const char *
        print_elem_type       ( )
        {
            return "double";
        }

}; // struct elem_type_traits <double>


/******************************
 ** struct: elem_type_traits **
 **     (spec: cxdb)         **
 ******************************/
template <>
struct elem_type_traits <cxdb>
{

    public:

        typedef cxdb elem_type;
        typedef double value_type;

        static inline
        const char *
        print_elem_type       ( )
        {
            return "cxdb";
        }

}; // struct elem_type_traits <cxdb>


# endif // __ELEM_TYPE_TRAITS_HPP__
