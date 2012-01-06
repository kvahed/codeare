AC_PREREQ([1.10])
AC_INIT([jrrs], [0.9], [k.vahedipour@fz-juelich.de])
AC_CONFIG_SRCDIR([src])

# m4 macros
AC_CONFIG_MACRO_DIR([config])
# where to put auxilliary files
AC_CONFIG_AUX_DIR([config])

AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE#([foreign])
LT_PREREQ([2.2])
LT_INIT([dlopen])

# Revision string
define([svnversion], esyscmd([sh -c "svnversion|tr -d '\n'"]))dnl
AC_DEFINE(SVN_REVISION, "svnversion", [SVN Revision])
# copyright notice
AC_COPYRIGHT([2010, 2011, Kaveh Vahedipour])

# Produce first version of idl derivates
esyscmd([sh -c "cd src; omniidl -bcxx RRSModule.idl; cd .."])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL

# Set programming language
AC_LANG([C++])

# Check SMP capabilities
AC_OPENMP

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_TIME
AC_HEADER_STDBOOL
AC_C_INLINE

AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# canonical host system type string
AC_CANONICAL_HOST

# Check for Python
m4_define([_AM_PYTHON_INTERPRETER_LIST], [python python3 python3.1])
AC_ARG_WITH([pythonpath],
AC_HELP_STRING([--with-pythonpath=PATH],
  [specify an absolute path to python executable when automatic version check (incorrectly) fails]),
  [python_path="$withval"], [python_path="auto"])
if test "$python_path" = auto; then
	AM_PATH_PYTHON([2.7])
else
	# Older versions of automake can't handle python3 well.  This is an
	# in-house workaround for them.
	PYTHON=$python_path
	AC_SUBST(PYTHON)
	PYTHON_PREFIX='${prefix}'
	AC_SUBST(PYTHON_PREFIX)
	PYTHON_EXEC_PREFIX='$(exec_prefix)'
	AC_SUBST(PYTHON_EXEC_PREFIX)
	PYTHON_VERSION=[`$PYTHON -c "import sys; sys.stdout.write(sys.version[:3])"`]
	if test `echo "$PYTHON_VERSION >= 2.7" | bc` != 1 ; then
		AC_MSG_ERROR(["Python version too old: $PYTHON_VERSION, need 2.7 or higher"])
	fi
	AC_SUBST(PYTHON_VERSION)
	PYTHON_PLATFORM=`$PYTHON -c "import sys; print(sys.platform)"`
	AC_SUBST(PYTHON_PLATFORM)
	pythondir='${prefix}/lib/python'$PYTHON_VERSION'/site-packages'
	AC_SUBST(pythondir)
	pkgpythondir='${pythondir}/'$PACKAGE
	AC_SUBST(pkgpythondir)
	pyexecdir='${exec_prefix}/lib/python'$PYTHON_VERSION'/site-packages'
	AC_SUBST(pyexecdir)
	pkgpyexecdir='${pyexecdir}/'$PACKAGE
	AC_SUBST(pkgpyexecdir)
fi

# Check for python development environments
if test -x ${PYTHON}-config; then
	PYTHON_INCLUDES=`${PYTHON}-config --includes`

	for flag in `${PYTHON}-config --ldflags`; do
		# add any '-L..." flags to PYTHON_LDFLAGS
		flag=`echo $flag | sed -ne 's/^\(\-L.*\)$/\1/p'`
		if test "X${flag}" != X; then
			PYTHON_LDFLAGS="$PYTHON_LDFLAGS ${flag}"
		fi
	done
	# on some platforms, ${PYTHON}-config --ldflags doesn't provide a -L
	# option while having the library under a non trivial directory.
	# as a workaround we try the "lib" sub directory under the common
	# prefix for this python.
	if test -z "${PYTHON_LDFLAGS}"; then
		PYTHON_LDFLAGS="-L`${PYTHON}-config --prefix`/lib"
	fi
else
	if test "X$PYTHON_INCLUDES" = X -o "X$PYTHON_LDFLAGS" = X; then
		AC_MSG_WARN([${PYTHON}-config does not exist or is not executable, so we could not detect python development environment.  Your system may require an additional package (e.g. "python3-dev").  Alternatively, if you are sure you have python headers and libraries, define PYTHON_INCLUDES and PYTHON_LDFLAGS and run this script.])
	fi
fi
AC_SUBST(PYTHON_INCLUDES)
AC_SUBST(PYTHON_LDFLAGS)
CXXFLAGS="${CXXFLAGS} ${PYTHON_INCLUDES}"

# enable "const" keyword
AC_C_CONST

# enable "restrict" keyword
AC_C_RESTRICT

# enable "inline" keyword
AC_C_INLINE

# Compiler vendor (Courtesy fftw)
AX_COMPILER_VENDOR

AC_CHECK_TYPES([ptrdiff_t])

#AC_CHECK_FUNCS([isnan isinf])
AC_DEFUN([AC_CXX_NAMESPACES],
	[AC_MSG_CHECKING([whether C++ namespaces are supported])
	 AC_LANG_PUSH([C++])
	 AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[namespace Outer { namespace Inner { int i = 0; }}]], [[using namespace Outer::Inner; return i;]])],[ac_cv_cxx_namespaces=yes],[ac_cv_cxx_namespaces=no])
     AC_LANG_POP
	 AC_MSG_RESULT($ac_cv_cxx_namespaces)
	 if test "$ac_cv_cxx_namespaces" = yes; then
	     AC_DEFINE(HAVE_NAMESPACES,,[define if the compiler implements namespaces])
     fi])
			      
AC_DEFUN([AC_CXX_ISNAN_IN_NAMESPACE_STD],
	[AC_REQUIRE([AC_CXX_NAMESPACES])
	 AC_MSG_CHECKING([whether isnan is in namespace std])
	 AC_LANG_PUSH([C++])
	 AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <cmath>
	   namespace lemur { int isnan(float x){ return std::isnan(x); } };]], [[
	   using namespace lemur; float x = 1.0; lemur::isnan(x); return 0;]])],[ac_cv_cxx_isnan_std=yes],[ac_cv_cxx_isnan_std=no])
     AC_LANG_POP
 	 AC_MSG_RESULT($ac_cv_cxx_isnan_std)
     if test "$ac_cv_cxx_isnan_std" = yes; then
	     AC_DEFINE(ISNAN_IN_NAMESPACE_STD,,[define if the compiler has isnan function in namespace std])
	 fi])

AC_DEFUN([AC_CXX_ISNAN_IN_NAMESPACE_GNU_CXX],
	[AC_REQUIRE([AC_CXX_NAMESPACES])
	 AC_MSG_CHECKING([whether isnan is in namespace __gnu_cxx])
	 AC_LANG_PUSH([C++])
	 AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <cmath>
	   namespace lemur { int isnan(float x){ return __gnu_cxx::isnan(x); } };]], [[
	   using namespace lemur; float x = 1.0; lemur::isnan(x); return 0;]])],[ac_cv_cxx_isnan_gnu_cxx=yes],[ac_cv_cxx_isnan_gnu_cxx=no])
     AC_LANG_POP
 	 AC_MSG_RESULT($ac_cv_cxx_isnan_gnu_cxx)
     if test "$ac_cv_cxx_isnan_gnu_cxx" = yes; then
	     AC_DEFINE(ISNAN_IN_NAMESPACE_GNU_CXX,,[define if the compiler has isnan function in namespace __gnu_cxx])
	 fi])

AC_CXX_ISNAN_IN_NAMESPACE_STD	       
AC_CXX_ISNAN_IN_NAMESPACE_GNU_CXX


AC_HEADER_MAJOR
# Checks for header files.
AC_CHECK_HEADERS([assert.h],         [], [AC_MSG_ERROR([cannot find assert.h, bailing out])])
#AC_CHECK_HEADERS([complex],          [], [AC_MSG_ERROR([cannot find complex, bailing out])])
AC_CHECK_HEADERS([ctype.h],          [], [AC_MSG_ERROR([cannot find ctype.h, bailing out])])
AC_CHECK_HEADERS([typeinfo],         [], [AC_MSG_ERROR([cannot find typeinfo, bailing out])])
AC_CHECK_HEADERS([fstream],          [], [AC_MSG_ERROR([cannot find fstream, bailing out])])
AC_CHECK_HEADERS([iostream],         [], [AC_MSG_ERROR([cannot find iostream, bailing out])])
AC_CHECK_HEADERS([limits.h],         [], [AC_MSG_ERROR([cannot find limits.h, bailing out])])
AC_CHECK_HEADERS([stdlib.h],         [], [AC_MSG_ERROR([cannot find stdio.h, bailing out])])
AC_CHECK_HEADERS([stdlib.h],         [], [AC_MSG_ERROR([cannot find stdlib.h, bailing out])])
AC_CHECK_HEADERS([string.h],         [], [AC_MSG_ERROR([cannot find string.h, bailing out])])
AC_CHECK_HEADERS([signal.h],         [], [AC_MSG_ERROR([cannot find signal.h, bailing out])])
AC_CHECK_HEADERS([omniORB4/CORBA.h], [], [AC_MSG_ERROR([cannot find omniORB4/CORBA.h, bailing out])])
AC_CHECK_HEADERS([nfft3.h],          [], [AC_MSG_RESULT([cannot find nfft3.h, NuFFT/CGSENSE unavailable!])])
AC_CHECK_HEADERS([nfft3util.h],      [], [AC_MSG_RESULT([cannot find nfft3util.h, NuFFT/CGSENSE unavailable!])])
AC_CHECK_HEADERS([H5Cpp.h],          [], [AC_MSG_RESULT([cannot find H5Cpp.h, hdf5 dump / read functionality unavailable!])])
AC_CHECK_HEADERS([cdf.h],            [], [AC_MSG_RESULT([cannot find cdf.h, cdf dump / read functionality unavailable!])])
AC_CHECK_HEADERS([nifti1_io.h],      [], [AC_MSG_RESULT([cannot find nifti1_io.h, NIFTI dump / read functionality unavailable!])])
AC_CHECK_HEADERS([gsl/gsl_wavelet2d.h], [HAVE_GSL=1], [AC_MSG_RESULT([cannot find gsl/gsl_wavelet2d.h, CompressedSensing unavailable!])])
AC_CHECK_HEADERS([sys/time.h c_asm.h intrinsics.h mach/mach_time.h])
AC_CHECK_HEADERS([fcntl.h],          [], [AC_MSG_ERROR([cannot find fcntl.h, bailing out])])
AC_CHECK_HEADERS([stddef.h],         [], [AC_MSG_ERROR([cannot find stddef.h, bailing out])])
AC_CHECK_HEADERS([cblas.h],          [], [AC_MSG_ERROR([cannot find cblas.h, bailing out])])

AC_CHECK_HEADERS(CL/cl.h, [AC_DEFINE(HAVE_OPENCL_HEADERS, 1, [Open CL headers])], [
    AC_CHECK_HEADERS(OpenCL/cl.h, [AC_DEFINE(HAVE_OPENCL_HEADERS, 1, [Open CL headers])], [
        echo "cannot file cl.h, GPU bloch simulator unavailable"
    ])
])

AC_MSG_CHECKING(whether host is Mac OS)
case $host in
    *darwin*)
        LDFLAGS="$LDFLAGS -framework OpenCL"
        AC_MSG_RESULT(yes)
    ;;
    *)
        AC_MSG_RESULT(no)
        AC_CHECK_LIB(OpenCL, clGetPlatformIDs, [],[ 
            echo "OpenCL library not found in your system"
            exit -1
        ])
esac

# Checks for C libraries.
AC_CHECK_LIB([omnithread],   [read],          [], [AC_MSG_ERROR([cannot find libomnithread, bailing out])])
AC_CHECK_LIB([omniDynamic4], [main],          [], [AC_MSG_ERROR([cannot find libomniDynamic, bailing out])])
AC_CHECK_LIB([dl],           [dlopen],        [], [AC_MSG_ERROR([cannot find libdl, bailing out])])
AC_CHECK_LIB([z],            [deflate],       [], [AC_MSG_ERROR([cannot find libz, bailing out])])
AC_CHECK_LIB([fftw3],        [fftw_execute],  [], [AC_MSG_ERROR([cannot find libfftw, bailing out])])
AC_CHECK_LIB([fftw3f],       [fftwf_execute], [], [AC_MSG_ERROR([cannot find libfftw3f, bailing out])])
AC_CHECK_LIB([cdf],          [CDFcreatezVar], [], [AC_MSG_RESULT([cannot find libcdf, cdf dump / read functionality unavailable!])])
AC_CHECK_LIB([gsl],          [gsl_wavelet_alloc], [], [AC_MSG_RESULT([cannot find libgsl, Compressed sensing reconstruction unavailable!])])
AC_CHECK_LIB([nfft3],        [nfft_precompute_psi], [HAVE_NFFT3=1], [AC_MSG_RESULT([cannot find libnfft3, NuFFT/CGSENSE unavailable!])])

AM_CONDITIONAL([USE_NFFT3],  [test "$HAVE_NFFT3" -eq 1])
AM_CONDITIONAL([USE_GSL],    [test "$HAVE_GSL" -eq 1])

# Time & Co for FFTW cycle
AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t is defined in <sys/time.h>])],,[
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
])
AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time])
dnl Cray UNICOS _rtc() (real-time clock) intrinsic
AC_MSG_CHECKING([for _rtc intrinsic])
rtc_ok=yes
AC_TRY_LINK([
#ifdef HAVE_INTRINSICS_H
#include <intrinsics.h>
#endif
], [_rtc()], [AC_DEFINE(HAVE__RTC,1,[Define if you have the UNICOS _rtc() intrinsic.])], [rtc_ok=no])
AC_MSG_RESULT($rtc_ok)


# Checks for C++ libraries with name spaces

# HDF5
AC_MSG_CHECKING([for getOrder in -lhdf5_cpp])
OLD_LIBS=$LIBS
LIBS="${LIBS} -lhdf5_cpp -lhdf5"

AC_TRY_LINK([#include <H5Cpp.h>],
   [
     H5::H5File file ("test.bin", H5F_ACC_RDONLY); 
	 H5::DataSet dataset = file.openDataSet("IntArray"); 
	 H5::IntType intype = dataset.getIntType(); 
	 H5std_string order_string; 
	 H5T_order_t order = intype.getOrder( order_string );
   ],
   [AC_MSG_RESULT([yes])]
   hdf5_found=yes,
   [LIBS=$OLD_LIBS]
   [AC_MSG_RESULT([no. Cannot find libhdf5_cpp, hdf5 dump/read not available.])])

if test -z $hdf5_found; then
    LIBS=$OLD_LIBS
fi


# NIFTI
AC_MSG_CHECKING([for nifti_image_read in -lniftiio])
OLD_LIBS=$LIBS
LIBS="${LIBS} -lniftiio -lznz"

AC_TRY_LINK([#include <nifti1_io.h>],
   [
	 const char* fname; nifti_image* ni = nifti_image_read (fname, 1);
   ],
   [AC_MSG_RESULT([yes])]
   niftiio_found=yes,
   [LIBS=$OLD_LIBS]
   [AC_MSG_RESULT([no. Cannot find libniftiio, NIFTI dump/read not available.])])

if test -z $niftiio_found; then
    LIBS=$OLD_LIBS
fi

# ITK
AH_TEMPLATE([HAVE_INSIGHT], [Define to 1 if you have the <itkImage.h> header file.])

AC_MSG_CHECKING([for InsightToolkit])
AC_ARG_WITH(itkroot,[  --with-itkroot       path to ITK install dir], [], [with_itkroot=no])
ITK_INSTALL_PREFIX=${with_itkroot}
ITK_INSTALL_PREFIX=`echo $ITK_INSTALL_PREFIX | sed -e 's|/|\\\/|g'`
itklib_path="${with_itkroot}/lib/InsightToolkit"
itkconfig_cmake="${itklib_path}/ITKConfig.cmake"
if test -a "$itkconfig_cmake"; then
	INC_INSIGHT_H="`cat $itkconfig_cmake | grep 'SET(ITK_INCLUDE_DIRS' | sed 's/SET(ITK_INCLUDE_DIRS \"/;/' | sed 's/\")//' | sed 's/;/ -I/g'`"
	INC_INSIGHT_H="`echo $INC_INSIGHT_H | sed -e "s|ITK_INSTALL_PREFIX|${ITK_INSTALL_PREFIX}|g" | sed -e 's|\${||g' | sed -e 's|}||g'`"
    insight_libs_dir="`cat $itkconfig_cmake | grep 'SET(ITK_LIBRARY_DIRS' | sed 's/SET(ITK_LIBRARY_DIRS \"/;/' | sed 's/\")//' | sed 's/;//g' | sed -e 's/\${ITK_INSTALL_PREFIX}/\/usr\/local/g'`"
	INSIGHT_LIB_DIR_FLAG="-L$insight_libs_dir"	
	INSIGHT_LIB_DIR_FLAG="`echo $INSIGHT_LIB_DIR_FLAG | sed -e "s|ITK_INSTALL_PREFIX|${ITK_INSTALL_PREFIX}|g" | sed -e 's|\${||g' | sed -e 's|}||g'`"
	INSIGHT_LIB_FLAG="-lITKFEM -lITKIO -lITKMetaIO -lITKAlgorithms -lITKStatistics -lITKBasicFilters -lITKCommon  -litkvnl -litkvnl_algo -litksys -lITKDICOMParser -lz"
    DEF_INSIGHT="#define HAVE_INSIGHT"
	AC_DEFINE(HAVE_INSIGHT, "1")
    HAVE_INSIGHT="yes"

    AC_MSG_RESULT([found.])
else
    AC_MSG_RESULT([not found.])
fi
         
AC_SUBST(ITK_INSTALL_PREFIX)
AC_SUBST(HAVE_INSIGHT)
AC_SUBST(INSIGHT_LIB_FLAG)
AC_SUBST(INSIGHT_LIB_DIR_FLAG)
AC_SUBST(INC_INSIGHT_H)

AC_MSG_CHECKING([for itk::Image in ITK])
OLD_LIBS=$LIBS
OLD_CXXFLAGS=$CXXFLAGS
LIBS="${LIBS} ${INSIGHT_LIB_DIR_FLAG} ${INSIGHT_LIB_FLAG}"
CXXFLAGS="${CXXFLAGS} ${INC_INSIGHT_H}"

AC_TRY_LINK([#include "itkImage.h"],
   [
	 typedef itk::Image< unsigned short, 3 > ImageType; ImageType::Pointer image = ImageType::New();
   ],
   [AC_MSG_RESULT([yes])]
   itk_found=yes,
   [LIBS=$OLD_LIBS]
   [AC_MSG_RESULT([no. Cannot find libITKCommon, ITK operations not available.])])

#if test -z $itk_found; then
    LIBS=$OLD_LIBS
	CXXFLAGS=$OLD_CXXFLAGS
#fi


# OMNIORB
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
    LIBS=$OLD_LIBS
fi

# MATLAB IO LIBRARIES
AC_ARG_WITH(matlabroot,[  --with-matlabroot       path to MATLAB install dir], 
        [], [with_matlabroot=no])

if test x"$with_matlabroot" != x"no" ; then

    MATLAB_ROOT="$with_matlabroot"
    MATLAB_LIB_PATH="${MATLAB_ROOT}/bin"
	CXXFLAGS="${CXXFLAGS} -I${MATLAB_ROOT}/extern/include"
	CPPFLAGS="${CPPFLAGS} -I${MATLAB_ROOT}/extern/include"

	AC_CHECK_HEADERS([mat.h],         [], [AC_MSG_ERROR([cannot find mat.h, bailing out])])

	AC_MSG_CHECKING([for MATFile in -lmat -lmx])

	OLD_LIBS=$LIBS

    case $host in
        *darwin*)
            MATLAB_LIB_PATH="${MATLAB_LIB_PATH}/maci64"
        ;;
		*x86_64*linux*)
		    MATLAB_LIB_PATH="${MATLAB_LIB_PATH}/glnxa64"
		;;
    esac

	LIBS="${LIBS} -L${MATLAB_LIB_PATH} -lmat -lmx -Wl,-rpath,${MATLAB_LIB_PATH}"

	AC_TRY_LINK([#include "mat.h"],
        [MATFile *pmat;],
   		[AC_MSG_RESULT([yes])]
   		matlab_found=yes,
   		[LIBS=$OLD_LIBS]
   		[AC_MSG_ERROR([no. cannot link to libmat, bailing out])])

fi

# BLAS and LAPACK
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

# If we're on configuring in OSX
if test -z "$BLAS"; then
case $host in
  *darwin*)
    BLAS="-DYA_BLAS -DYA_LAPACK -DYA_BLASMULT -lblas -lcblas -llapack -framework vecLib"
    AUTOAPPLE="yes"
    AC_MSG_RESULT([Mac OSX: Assuming vecLib is available!])
    ;;
esac
fi

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

AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([getpagesize])
AC_CHECK_FUNCS([memmove])
AC_CHECK_FUNCS([memset])
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([setenv])
AC_CHECK_FUNCS([sqrt])
AC_CHECK_FUNCS([strchr])
AC_CHECK_FUNCS([strstr])
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MMAP


# Generate files
AC_CONFIG_FILES([src/Loader.hpp Makefile doc/Makefile src/Makefile src/modules/Makefile src/tinyxml/Makefile])
AC_OUTPUT

