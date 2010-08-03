AC_PREREQ([1.10])
AC_INIT([jrrs], [0.4], [k.vahedipour@fz-juelich.cd])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE([foreign])
LT_INIT([dlopen])
# Revision string
define([svnversion], esyscmd([sh -c "svnversion|tr -d '\n'"]))dnl
AC_DEFINE(SVN_REVISION, "svnversion", [SVN Revision])

# Produce first version of idl derivates
esyscmd([sh -c "cd src; omniidl -bcxx RRSModule.idl; cd .."])

# Checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL

# Set programming language
AC_LANG(C++)

# Checks for header files.
AC_CHECK_HEADERS([complex],          [], [AC_MSG_ERROR([cannot find complex, bailing out])])
AC_CHECK_HEADERS([assert.h],         [], [AC_MSG_ERROR([cannot find assert.h, bailing out])])
AC_CHECK_HEADERS([stdlib.h],         [], [AC_MSG_ERROR([cannot find cstdlib, bailing out])])
AC_CHECK_HEADERS([string.h],         [], [AC_MSG_ERROR([cannot find ctime, bailing out])])
AC_CHECK_HEADERS([iostream],         [], [AC_MSG_ERROR([cannot find iostream, bailing out])])
AC_CHECK_HEADERS([fstream],          [], [AC_MSG_ERROR([cannot find fstream, bailing out])])
AC_CHECK_HEADERS([signal.h],         [], [AC_MSG_ERROR([cannot find signal.h, bailing out])])
AC_CHECK_HEADERS([omniORB4/CORBA.h], [], [AC_MSG_ERROR([cannot find omniORB4/CORBA.h, bailing out])])

# Checks for C libraries.
AC_CHECK_LIB([omnithread],   [read],   [], [AC_MSG_ERROR([cannot find libomnithread, bailing out])])
AC_CHECK_LIB([omniDynamic4], [main],   [], [AC_MSG_ERROR([cannot find libomniDynamic, bailing out])])
AC_CHECK_LIB([dl],           [dlopen], [], [AC_MSG_ERROR([cannot find libdl, bailing out])])

# Checks for C++ libraries with name spaces
AC_MSG_CHECKING([for CORBA::ORB_init(int, char**) in -lomniORB4])
OLD_LIBS=$LIBS
LIBS="${LIBS} -lomniORB4"
cd src
omniidl -bcxx RRSModule.idl
cd ..
AC_TRY_LINK([#include "src/RRSModule.hh"],
   [int argc = 0; char** argv = 0; CORBA::ORB_var orb = CORBA::ORB_init(argc, argv);],
   [AC_MSG_RESULT([yes])]
   omniORB4_found=yes,
   [LIBS=$OLD_LIBS]
   [AC_MSG_ERROR([no. cannot find libomniORB4, bailing out])])
if test -z $omniORB4_found; then
    LIBS="${LIBS} -lomniORB4"
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CONFIG_FILES([Makefile doc/Makefile src/Makefile src/modules/Makefile src/Loader.h])
AC_OUTPUT

