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
AC_LANG([C++])
AC_OPENMP

# Checks for header files.
AC_CHECK_HEADERS([assert.h],         [], [AC_MSG_ERROR([cannot find assert.h, bailing out])])
AC_CHECK_HEADERS([complex],          [], [AC_MSG_ERROR([cannot find complex, bailing out])])
AC_CHECK_HEADERS([ctype.h],          [], [AC_MSG_ERROR([cannot find ctype, bailing out])])
AC_CHECK_HEADERS([fstream],          [], [AC_MSG_ERROR([cannot find fstream, bailing out])])
AC_CHECK_HEADERS([iostream],         [], [AC_MSG_ERROR([cannot find iostream, bailing out])])
AC_CHECK_HEADERS([stdlib.h],         [], [AC_MSG_ERROR([cannot find cstdlib, bailing out])])
AC_CHECK_HEADERS([string.h],         [], [AC_MSG_ERROR([cannot find ctime, bailing out])])
AC_CHECK_HEADERS([signal.h],         [], [AC_MSG_ERROR([cannot find signal.h, bailing out])])
AC_CHECK_HEADERS([omniORB4/CORBA.h], [], [AC_MSG_ERROR([cannot find omniORB4/CORBA.h, bailing out])])
AC_CHECK_HEADERS([noncart/nufft.h],  [], [AC_MSG_RESULT([cannot find noncart/nufft.h, NuFFT/CGSENSE unavailable!])])

#AC_CHECK_HEADERS([cblas.h],          [], [AC_MSG_RESULT([cannot find cblas.h, blas functions not available])])

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

AC_MSG_CHECKING([for noncart::nufft() in -lnoncart])
OLD_LIBS=$LIBS
LIBS="${LIBS} -lnoncart"
AC_TRY_LINK([#include "noncart/nufft.h"],
   [int argc = 0; noncart::nufft nufft;],
   [AC_MSG_RESULT([yes])]
   noncart_found=yes,
   [LIBS=$OLD_LIBS]
   [AC_MSG_RESULT([no. cannot find libnoncart, NuFFT/CGSENSE unavailable ])])
AM_CONDITIONAL([NONCART], [test x$noncart_found = xyes])

AC_ARG_WITH(blas,[  --with-blas             which blas library to use],[
case $withval in
  yes)
    BLAS='-lblas'
    ;;
  no)
    BLAS=''
    ;;
  *)
    BLAS=$withval
    ;;
esac
])

AC_ARG_WITH(lapack,[  --with-lapack           which lapack library to use],[
case $withval in
  yes)
    LAPACK='-llapack'
    ;;
  no)
    LAPACK=''
    ;;
  *)
    LAPACK=" $withval"
    ;;
esac
])

DGEMM_SYMBOL=cgemm_
DGEEV_SYMBOL=cgeev_

SAVE_LIBS=$LIBS
LIBS="$LIBS $CLIBS $FLIBS"

# If we're using the Intel compilers, we prefer MKL as well                                                                                 
if test -z "$BLAS"; then
  if test $CC = icc -o $CXX = icpc; then
    AC_CHECK_LIB(mkl, $DGEMM_SYMBOL, BLAS="-lmkl")
  fi
fi

# Otherwise, run through a list of possibilities                                                                                            
if test -z "$BLAS"; then
  AC_CHECK_LIB(essl, $DGEMM_SYMBOL, BLAS="-lessl -lpthread")
fi

if test -z "$BLAS"; then
  AC_CHECK_LIB(scs, $DGEMM_SYMBOL, BLAS="-lscs")
fi

if test -z "$BLAS"; then
  SAVE_LIBS2=$LIBS
  LIBS="-lpthread $LIBS"
  AC_CHECK_LIB(goto, $DGEMM_SYMBOL, BLAS="-lgoto -lpthread")
  LIBS=$SAVE_LIBS2
fi

if test -z "$BLAS"; then
  SAVE_LIBS2=$LIBS
  LIBS="-latlas $LIBS"
  AC_CHECK_LIB(f77blas, $DGEMM_SYMBOL, BLAS="-lcblas -latlas")
  LIBS=$SAVE_LIBS2
fi

if test -z "$BLAS"; then
  AC_CHECK_LIB(mkl, $DGEMM_SYMBOL, BLAS="-lmkl")
fi

if test -z "$BLAS"; then
case $target_vendor in
  apple)
    BLAS="-altivec -framework vecLib"
    AUTOAPPLE="yes"
    AC_MSG_RESULT([No BLAS specified and vendor is Apple: Assuming vecLib is available!])
    ;;
esac
fi

if test -z "$BLAS"; then
  AC_CHECK_LIB(blas, $DGEMM_SYMBOL, BLAS="-lblas")
fi

if test -z "$BLAS"; then
  AC_MSG_ERROR("Did not find a BLAS library.")
fi
AC_SUBST(BLAS)
#AC_DEFINE(HAVE_BLAS)
AC_MSG_RESULT([BLAS library is $BLAS])
LIBS=$SAVE_LIBS

# Similar checks for LAPACK, except that sometimes BLAS includes LAPACK                                                                     
# Also note that LAPACK often (always?) depends on BLAS                                                                                     
SAVE_LIBS=$LIBS
LIBS="$BLAS $LIBS $CLIBS $FLIBS"
if test -z "$LAPACK"; then
  if test ! -z "$AUTOAPPLE"; then
    AC_MSG_RESULT([Assuming vecLib for LAPACK])
    LAPACK="$BLAS"
    LAPACK_IN_BLAS=yes
  else
    AC_CHECK_FUNC($DGEEV_SYMBOL,LAPACK_IN_BLAS=yes)
  fi
fi
if test -z "$LAPACK_IN_BLAS"; then
  if test -z "$LAPACK"; then
    if test $CC = icc -o $CXX = icpc; then
      # Intel's libmkl_lapack64 contains the double-precision routines                                                                      
      AC_CHECK_LIB(mkl_lapack64,$DGEEV_SYMBOL,LAPACK="-lmkl_lapack64")
    fi
  fi
  if test -z "$LAPACK"; then
    AC_CHECK_LIB(lapack,$DGEEV_SYMBOL,LAPACK="-llapack")
  fi
  if test -z "$LAPACK"; then
    AC_CHECK_LIB(lapack-3,$DGEEV_SYMBOL,LAPACK="-llapack-3")
  fi
  if test -z "$LAPACK"; then
    AC_MSG_ERROR("Did not find a LAPACK library.");
  fi
  AC_MSG_RESULT([LAPACK library is $LAPACK])
else
  AC_MSG_RESULT([LAPACK library is included in BLAS+LIBS])
fi
AC_SUBST(LAPACK)
#AC_DEFINE(HAVE_LAPACK)
LIBS=$SAVE_LIBS


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CONFIG_FILES([Makefile doc/Makefile src/Makefile src/modules/Makefile src/Loader.h])
AC_OUTPUT

